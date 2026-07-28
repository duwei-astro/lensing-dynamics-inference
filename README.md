# lensing-dynamics-inference

Code and data for hierarchical inference of the post-Newtonian parameter
$\gamma_{\rm PPN}$ and stellar orbital anisotropy from galaxy-scale strong
lenses.

## Contents

```text
lensing-dynamics-inference/
├── README.md
├── LICENSE
├── requirements.txt
├── notebooks/
│   ├── run_mcmc_bpl_fiducial.ipynb
│   ├── *.py
│   ├── figures/
│   │   └── plot_*.ipynb
│   └── supplementary_figures/
│       └── plot_*.ipynb
└── data/
    ├── xbeta_likelihood_maps/
    │   └── *.fits
    ├── Fig3_beta_posterior_bands.csv
    ├── nsource_flag.dat
    ├── obs125_krms_vbias_kenv.dat
    └── Supplementary_table3_lens_sample.dat
```

The repository contains:

* an example notebook implementing the fiducial hierarchical inference and its supporting Python modules;
* notebooks used to reproduce the main-text and Supplementary figures;
* the $X$ – $\beta$ likelihood maps for all 125 modelled lenses;
* small supporting data used in the analysis and figure production.

The likelihood files include the four outliers
`BELLSJ0151+0049`, `BELLSJ1159-0007`, `BELLSJ1215+0047`, and
`BELLSJ1352+3216`. These systems were modelled and have corresponding
$X$ – $\beta$ maps, but were excluded from the final hierarchical
Bayesian inference, which was based on 121 lenses.

The supporting data include:

* `obs125_krms_vbias_kenv.dat`: input data used in the hierarchical inference, including rms and environmental convergence;
* `Fig3_beta_posterior_bands.csv`: posterior bands for the inferred redshift evolution of stellar orbital anisotropy shown in Fig. 3;
* `Supplementary_table3_lens_sample.dat`: source data for Supplementary Table 3, containing the observed properties and modelling inputs;
* `nsource_flag.dat`: number of background-source components adopted in the lens modelling for each system.

The hierarchical-inference chains used for the reported results and the
lens-modelling MCMC chains for the lenses shown in Supplementary Fig. 1 are
archived on Zenodo at https://doi.org/10.5281/zenodo.21188730.

## System requirements

The code requires Python 3.11 or later and has been tested with:

* Python 3.13.9;
* Rocky Linux 9.3 (Blue Onyx), x86_64;
* the Python dependencies and package versions listed in
  `requirements.txt`.

No non-standard hardware is required.

## Installation

Clone the repository and install the required Python packages:

```bash
git clone https://github.com/duwei-astro/lensing-dynamics-inference.git
cd lensing-dynamics-inference
python -m pip install -r requirements.txt
```

Installation typically takes a few minutes.

## Running the example notebook

Open

```text
notebooks/run_mcmc_bpl_fiducial.ipynb
```
The notebook implements the fiducial hierarchical inference used in the manuscript.

For a quick test of the workflow, users may reduce `nwalkers`, `nsteps`, and
the PSO parameter `n_particles`. The numbers of parallel processes specified
in `Pool(...)` and by `n_processes` should be adjusted according to the
available computational resources.

For a reduced demonstration with

```text
n_particles = 10
n_processes = 20
nwalkers = 20
nsteps = 30
Pool(20)
```

Run all notebook cells in order. On the tested system, the reduced demonstration takes approximately 1 minute using 20 CPU cores. This configuration is intended only to demonstrate the workflow and does not provide converged posterior constraints. A successful run produces diagnostic output and an MCMC chain file, which is saved in the same directory as the notebook

The fiducial external-convergence-prior setting is `lambda_val = 1.0` and users may 
change `lambda_val` to explore the sensitivity of the inference to the weight assigned to the external-convergence priors. 
To test the sensitivity of the results to the adopted cosmology, users can vary `Omatter` in `cosmo_kb_func.data_collection`, for example by setting `Omatter=0.35`.

## Reproducing the figures

The figure-production notebooks are located in

```text
notebooks/figures/
notebooks/supplementary_figures/
```

They use the supporting data included in this repository and the production
MCMC chains in the Zenodo archive linked above.

The notebooks currently read the chain files from `../../../NA_zenodo/`.
Download the required files to this directory, or update the input paths
in the notebooks to their local locations.

Run all cells in the relevant notebook to reproduce the corresponding
main-text or Supplementary figure. The generated figures are saved in the
same directory as the notebook.

## Applying the inference to another lens sample

The framework provides a scalable way to extract population-level
information from strong-lens samples. To analyse another lens sample,
replace the likelihood maps and lens metadata with files in the same
formats as the provided examples. If their filenames or locations differ,
update the corresponding input paths in the example notebook.

## License

The source code in this repository is licensed under the BSD 3-Clause
“New” or “Revised” License. See `LICENSE` for details.

<!--
## Citation

When using this code or the accompanying data products, please cite:

> Wei Du et al., “[Paper title]”, [journal and DOI].

and the associated Zenodo archive:

> [Zenodo citation and DOI]

## Contact

For questions about the code or data, please contact:

Wei Du  
[Email address]
-->


