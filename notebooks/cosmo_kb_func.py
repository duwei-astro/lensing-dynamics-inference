import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, z_at_value
from astropy import units as u
from astropy import constants as const
from scipy.interpolate import RegularGridInterpolator

def load_lensing_data(manifest_path, fits_dir):
    """  
    Read a manifest file containing lens system information and load corresponding FITS data files.  

    This version can:  
    1. Automatically filter out comment lines starting with '#'.  
    2. Only read the first six columns of the file.  
    3. Robustly read data from the main HDU or the first extension in the FITS file.  

    Args:  
        manifest_path (str): Path to the '.dat' or '.txt' manifest file.  
        fits_dir (str): Path to the directory containing 125 FITS files.  

    Returns:  
        dict: A dictionary containing all loaded data, or None if an error occurs.  
    """ 
    # --- 1. Read the manifest file ---
    try:
        # Define column names for the six needed columns
        col_names = ['lens_name', 'lens_redshift', 'source_redshift', 'kext_rms', 'vbias', 'kenv'] #for obs125_krms_vbias_kenv.dat

        # Use pandas to read, with comment='#' to ignore lines starting with '#'
        manifest_df = pd.read_csv(
            manifest_path,
            sep=r'\s+',
            header=None,
            names=col_names,
            usecols=[0, 1, 2, 3, 4, 5], 
            comment='#'           # <--- Ignore lines starting with '#'
        )

        print(f"Successfully loaded manifest file: {manifest_path}")
        print("Automatically ignored comment lines starting with '#', and loaded only the first 6 columns.")
        print(f"Found {len(manifest_df)} lens systems.")

    except FileNotFoundError:
        print(f"Error: Manifest file not found at {manifest_path}")
        return None
    except Exception as e:
        print(f"Error occurred while reading manifest file: {e}")
        return None

    # --- 2. Build FITS file paths ---
    fits_filepaths = [os.path.join(fits_dir, f"{row.lens_name}.fits") for index, row in manifest_df.iterrows()]

    # --- 3. Read the 125 FITS files ---
    all_pdfs = []
    x_points = []
    y_points = []
    print(f"\nStarting to load {len(fits_filepaths)} FITS files...")

    for i, filepath in enumerate(fits_filepaths):
        try:
            with fits.open(filepath) as hdul:
                data = None
                if hdul[0].data is not None:
                    data = hdul[0].data
                    xdata = hdul[1].data
                    ydata = hdul[2].data
                else:
                    raise IOError("No readable data found in FITS file.")

                if data.shape != (128, 128):
                    print(f"Warning: File {os.path.basename(filepath)} has shape {data.shape}, expected (128, 128).")

                all_pdfs.append(data)
                x_points.append(xdata)
                y_points.append(ydata)

        except FileNotFoundError:
            print(f"Error: FITS file not found: {filepath}")
            all_pdfs.append(None)
            continue
        except Exception as e:
            print(f"Error while reading FITS file {filepath}: {e}")
            all_pdfs.append(None)
            continue

    successful_loads = sum(p is not None for p in all_pdfs)
    print(f"Successfully loaded {successful_loads} / {len(fits_filepaths)} PDFs.")

    # --- 4. Store all data in a dictionary and return ---
    lensing_data = {
        'manifest_df': manifest_df,
        'pdf_data': all_pdfs,
        'x_data': x_points,
        'y_data': y_points
    }

    return lensing_data


def get_cosmology_bins(cosmology, z_final: float, num_bins: int, verbose=False):
    """
    Compute redshift range, and the proper radial distance in each redshift bin.

    Returns:
        DataFrame: Each row provides z_start, z_end, z_median, and the proper travel distance for each bin.
    """
    total_lookback_time = cosmology.lookback_time(z_final)
    total_comoving_distance = cosmology.comoving_distance(z_final)

    if verbose:
        print(f"Cosmology model in use: {cosmology.name}")
        print(f"Custom parameters: H0={cosmology.H0}, Om0={cosmology.Om0}")
        print(f"Total lookback time to z={z_final} is {total_lookback_time:.3f}")
        print(f"Total comoving distance to z={z_final} is {total_comoving_distance:.3f}\n")
        print(f"The number of redshift bins is {num_bins:.1f}")

    lookback_time_bin_edges = np.linspace(0 * u.Gyr, total_lookback_time, num_bins + 1)

    redshift_bin_edges = []
    for t in lookback_time_bin_edges:
        if t.value == 0:
            redshift_bin_edges.append(0.0)
        else:
            z_quantity = z_at_value(cosmology.lookback_time, t)
            redshift_bin_edges.append(z_quantity.value)

    time_start = lookback_time_bin_edges[0:-1]
    time_end = lookback_time_bin_edges[1:]
    z_start = np.array(redshift_bin_edges[0:-1])
    z_end = np.array(redshift_bin_edges[1:])
    z_median = 0.5*(z_start+z_end)

    dist_start = cosmology.comoving_distance(z_start)
    dist_end = cosmology.comoving_distance(z_end)

    # This difference gives the actual travel path length within each bin
    light_path_distance = (dist_end - dist_start)/(1+z_median)
    proper_distance = (time_end - time_start)*306.60139371446905 #Mpc

    if verbose:
        print(f"The proper distances are {proper_distance[0:5].value} Mpc")

    # Combine all data as columns in a DataFrame
    redshift_df = pd.DataFrame({
        'z_start': z_start,
        'z_end': z_end,
        'z_med': z_median,
        'delta_distance': proper_distance.value  # unit number (Mpc)
    })

    return redshift_df

def get_weighting(cosmo,z_d,z_s,z_i,delta_dp):
    # get the strong lensing weighting of the mass bins at redshifts z_i
    if not hasattr(delta_dp, "unit") or not delta_dp.unit.is_equivalent(u.Mpc):
        raise ValueError("delta_dp must have units of u.Mpc, for example delta_dp = 1.2 * u.Mpc")
    H0 = cosmo.H0 #km/s/Mpc
    c = const.c.to(u.km / u.s) #km/s
    Om0 = cosmo.Om0
    D_s = cosmo.angular_diameter_distance(z_s) #Mpc
    D_i = cosmo.angular_diameter_distance(z_i) #Mpc
    z_i = np.minimum(z_s,z_i) # setting the weight to zero for bins with redshift greater than z_s
    D_is = cosmo.angular_diameter_distance_z1z2(z_i, z_s) #Mpc

    z_1 = np.minimum(z_d,z_i)
    z_2 = np.maximum(z_d,z_i)
    D_1s = cosmo.angular_diameter_distance_z1z2(z_1, z_s) #Mpc
    D_2 = cosmo.angular_diameter_distance(z_2) #Mpc
    D_12 = cosmo.angular_diameter_distance_z1z2(z_1, z_2) #Mpc

    bet = D_12/D_2*D_s/D_1s
    bet[z_i<z_d] = 0.0 #For EPL model; Johnson et al. 2025, Foreground_biases_in_SL
    w = (1-bet)*Om0*3*H0**2/2/c**2*D_is/D_s*D_i*(1+z_i)**3*delta_dp

    return w


def get_weight_list(cosmo, z_d, z_s, z_i, delta_dp):
    """
    Compute the lensing weight list for multiple lens systems.

    The function loops over each (z_d, z_s) pair and computes the weights at redshift slices z_i.

    Args:
        cosmo: astropy.cosmology object.
        z_d (array_like): 1D array of lens redshifts.
        z_s (array_like): 1D array of source redshifts. Must be the same length as z_d.
        z_i (array_like): Mid-redshift value of each bin (shared).
        delta_dp (array_like or astropy.Quantity): Radial proper distance thickness or light-travel-distance for each bin.

    Returns:
        list
    """
    # --- 1. Input validation ---
    if len(z_d) != len(z_s):
        raise ValueError("Input error: z_d and z_s must have the same number of entries.")

    # Ensure delta_dp is an astropy Quantity with units
    if not hasattr(delta_dp, "unit") or not delta_dp.unit.is_equivalent(u.Mpc):
        raise ValueError("delta_dp must have units of u.Mpc, for example delta_dp = 1.2 * u.Mpc")

    # --- 2. Main calculation --- Loop over all systems
    weights_list = []  
    for i in range(len(z_d)):
        wgt = get_weighting(cosmo,z_d[i],z_s[i],z_i,delta_dp)
        weights_list.append(wgt.values)

    return weights_list


# Combine data + create interpolators
def data_collection(manifest_path, fits_dir, num_systems=125, selected_indices=None,Omatter=0.315,Hzero=67.4):  
    """  
    Combine lensing system data and create interpolators.  
    
    Args:  
        manifest_path (str): Path to the manifest file.  
        fits_dir (str): Directory containing FITS files.  
        num_systems (int): Total number of systems to load (if selected_indices is None).  
        selected_indices (list or array-like, optional):   
            List or array of indices to select systems.  
            If None, select the first num_systems systems.  
            Default is None.  
            
    Returns:  
        dict: Dictionary containing processed data.  
    """ 
    print("--- Loading data ---")  
    # Read the redshifts and P(X, \beta) info for all 125 systems  
    lensing_data = load_lensing_data(manifest_path, fits_dir)  
    
    # If selected_indices is given, filter the data accordingly 
    if selected_indices is not None:  
        print(f"--- Select data by index: {selected_indices} ---")  
        # Ensure indices is array-like for indexing  
        indices = np.array(selected_indices)  
        
        # Use .iloc to select data from pandas Series by integer positions   
        z_d = lensing_data["manifest_df"]["lens_redshift"].iloc[indices].reset_index(drop=True)  
        z_s = lensing_data["manifest_df"]["source_redshift"].iloc[indices].reset_index(drop=True)  
        kext_rms = lensing_data["manifest_df"]["kext_rms"].iloc[indices].reset_index(drop=True)  
        vbias = lensing_data["manifest_df"]["vbias"].iloc[indices].reset_index(drop=True)  
        kenv = lensing_data["manifest_df"]["kenv"].iloc[indices].reset_index(drop=True)

        # For list-like data, use list comprehensions for filtering  
        pdf_data = [lensing_data["pdf_data"][i] for i in indices]  
        x_data = [lensing_data["x_data"][i] for i in indices]  
        y_data = [lensing_data["y_data"][i] for i in indices]  
        
        # Update num_systems to the length of selected systems  
        num_systems = len(indices)  
    else:  
        # If selected_indices is not given, default to selecting the first num_systems  
        print(f"--- Selecting first {num_systems} entries ---")   
        z_d = lensing_data["manifest_df"]["lens_redshift"][:num_systems]  
        z_s = lensing_data["manifest_df"]["source_redshift"][:num_systems]  
        kext_rms = lensing_data["manifest_df"]["kext_rms"][:num_systems]  
        vbias = lensing_data["manifest_df"]["vbias"][:num_systems] 
        kenv = lensing_data["manifest_df"]["kenv"][:num_systems]  
        pdf_data = lensing_data["pdf_data"][:num_systems]  
        x_data = lensing_data["x_data"][:num_systems]  
        y_data = lensing_data["y_data"][:num_systems]  

    cosmology = FlatLambdaCDM(H0=Hzero, Om0=Omatter, name='flat_lambda_CDM')  
    final_redshift = 1.52  
    number_of_bins = 50  

    cosmology_bins_df = get_cosmology_bins(cosmology, final_redshift, number_of_bins,verbose=True)  
    z_i = cosmology_bins_df["z_med"]  
    delta_Dp = cosmology_bins_df["delta_distance"].values*u.Mpc  
    
    # Generate weights (this will be based on the filtered z_d and z_s)  
    weights = get_weight_list(cosmology, z_d, z_s, z_i, delta_Dp)  
    
    # Create the interpolators for likelihood calculation  
    pdf_interpolators = []  
    # The loop count is determined by the updated num_systems   
    for i in range(num_systems):      
        pdf_grid = pdf_data[i]  
        kappa_range = x_data[i]  
        beta_range = y_data[i]  
        
        # Create and store the interpolator  
        interpolator = RegularGridInterpolator((beta_range, kappa_range), pdf_grid, bounds_error=False, fill_value=0)  
        pdf_interpolators.append(interpolator)  

    # Bundle all filtered data and return  
    dummy_data = {  
        "z_d": z_d,  
        "z_s": z_s,  
        "z_i": z_i,  
        "weights": weights,  
        "kext_rms": kext_rms,  
        "vbias": vbias,  
        "kenv": kenv,
        "pdf_interpolators": pdf_interpolators  
    }  
    print(f"--- Data packaging complete ({num_systems} groups in total) ---\n")  
    return dummy_data  


# Function to calculate the Highest Density Interval (HDI)
def hdi_interval(samples, credible_interval=0.68):
    """
    Compute the Highest Density Interval (HDI).
    
    Parameters:
    -----------
    samples : array-like
        Sample data
    credible_interval : float
        Credible interval (0-1)
    
    Returns:
    --------
    lower, upper : float
        Lower and upper bounds of the HDI
    """
    # Ensure samples is a 1D array
    samples = np.asarray(samples).flatten()
    samples = np.sort(samples)
    n = len(samples)
    
    # Compute the number of samples to include
    n_to_include = int(np.ceil(credible_interval * n))
    
    # Find the interval with minimum width
    interval_width = []
    interval_indices = []
    
    for i in range(n - n_to_include + 1):
        width = samples[i + n_to_include - 1] - samples[i]
        interval_width.append(width)
        interval_indices.append((i, i + n_to_include - 1))
    
    # Select the interval with minimum width
    min_idx = np.argmin(interval_width)
    lower_idx, upper_idx = interval_indices[min_idx]
    
    return samples[lower_idx], samples[upper_idx]

def compute_hdi_bands(y_pred_samples, credible_intervals=[0.68, 0.95]):
    """
    Compute HDI credible intervals for each predicted point.
    
    Parameters:
    -----------
    y_pred_samples : array, shape (n_samples, n_pred_points)
        Prediction samples
    credible_intervals : list
        List of credibility levels
    
    Returns:
    --------
    bands : dict
        Contains the lower and upper bounds for each credibility level
    """
    n_pred_points = y_pred_samples.shape[1]
    bands = {}
    
    for ci in credible_intervals:
        lower_bounds = np.zeros(n_pred_points)
        upper_bounds = np.zeros(n_pred_points)
        
        for i in range(n_pred_points):
            lower, upper = hdi_interval(y_pred_samples[:, i], ci)
            lower_bounds[i] = lower
            upper_bounds[i] = upper
        
        bands[f'{int(ci*100)}%'] = (lower_bounds, upper_bounds)
    
    return bands


def get_mcmc_title(mcsamples,percent=[16, 50, 84],symmetry_threshold = 1e-3,weight=False, fmt=None):
    """
    symmetry_threshold: Threshold for symmetry check
        If the absolute difference between upper and lower errors is below this threshold,
        treat the errors as effectively symmetric
    fmt: Format specifier(s) for each parameter.
        - None: use ".2f" for all parameters (default)
        - str: apply the same format to all parameters (e.g. ".3f")
        - list of str: per-parameter format (e.g. [".3f", ".2f", ".2f"]);
          if the list is shorter than ndim, ".2f" is used for remaining parameters
    """
    stats = mcsamples.getMargeStats()
    names = mcsamples.paramNames.list()
    ndim = len(names)
    txt_title = []
    perc_all = []
    
    for i in range(ndim):
        param = stats.parWithName(names[i]) 
        if weight:
            perc = np.percentile(mcsamples.samples[:,i],percent,weights=mcsamples.weights,method='inverted_cdf')
        else:
            perc = np.percentile(mcsamples.samples[:,i],percent)
        
        perc_all.append(perc.copy())
        q = np.diff(perc)  

        # Resolve fmt_specifier for this parameter
        if fmt is None:
            fmt_specifier = ".2f"
        elif isinstance(fmt, str):
            fmt_specifier = fmt
        else:  # list
            fmt_specifier = fmt[i] if i < len(fmt) else ".2f"

        # Use the selected format string
        central_value_str = f"{perc[1]:{fmt_specifier}}"

        # Check if errors are (approximately) symmetric
        if abs(q[0] - q[1]) < symmetry_threshold:
            # Symmetric case: use ±
            avg_error = (q[0] + q[1]) / 2.0
            error_str = f"{avg_error:{fmt_specifier}}"
            txt = f"${param.label} = {central_value_str} \\pm {error_str}$"
        else:
            # Asymmetric case: use sub/superscripts
            lower_err_str = f"{q[0]:{fmt_specifier}}"
            upper_err_str = f"{q[1]:{fmt_specifier}}"
            txt = f"${param.label} = {central_value_str}_{{ -{lower_err_str} }}^{{ +{upper_err_str} }}$"
            
        txt_title.append(txt)

    perc_all = np.asarray(perc_all)  # shape: (ndim, len(percent))
    return perc_all, txt_title

def merge_lens_data(d1, d2):
    combined = {}
    
    # Iterate through all keys
    for key in d1.keys():
        val1 = d1[key]
        val2 = d2[key]
        
        # 1. If it is a Pandas Series (z_d is this type)
        if isinstance(val1, pd.Series):
            # ignore_index=True is important, otherwise the indices will overlap (become 0..120, 0..120)
            combined[key] = pd.concat([val1, val2], ignore_index=True)
            
        # 2. If it is a regular Python list (corresponding to pdf_interpolators)
        elif isinstance(val1, list):
            combined[key] = val1 + val2
            
        # 3. If it is a NumPy array
        elif isinstance(val1, np.ndarray):
            combined[key] = np.concatenate([val1, val2])
            
        # 4. Other cases (e.g. a single value), either raise an error or keep the original
        else:
            print(f"Warning: Key '{key}' type {type(val1)} not handled. Keeping data1.")
            combined[key] = val1
            
    return combined

def inspect_data_struct(data):
    #Check the type, length/shape, and other information of all fields in the data dictionary
    import numpy as np
    import pandas as pd

    for key, val in data.items():
        # Type name
        t = type(val)

        # Try to get the length
        try:
            l = len(val)
        except TypeError:
            l = "No len()"

        # Try to get the shape
        shape = getattr(val, "shape", None)

        print(f"{key}:")
        print(f"  type   = {t}")
        print(f"  len    = {l}")
        print(f"  shape  = {shape}")
        print("-" * 40)

