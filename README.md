# lensing-dynamics-inference

Code and data for hierarchical inference of the post-Newtonian parameter
$\gamma_{\rm PPN}$ and stellar orbital anisotropy from galaxy-scale strong
lenses.

## Contents

```text
lensing-dynamics-inference/
├── README.md
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

## Running the example notebook

Download the repository and open

```text
notebooks/run_mcmc_bpl_fiducial.ipynb
```
The notebook implements the fiducial hierarchical inference used in the manuscript.

For a quick test of the workflow, users may reduce `nwalkers`, `nsteps`, and
the PSO parameter `n_particles`. The numbers of parallel processes specified
in `Pool(80)` and by `n_processes` should be adjusted according to the
available computational resources.

The fiducial external-convergence-prior setting is `lambda_val = 1.0` and users may 
change `lambda_val` to explore the sensitivity of the inference to the weight assigned to the external-convergence priors. 
To test the sensitivity of the results to the adopted cosmology, users can vary `Omatter` in `cosmo_kb_func.data_collection`, for example by setting `Omatter=0.35`.

The figure-production notebooks are located in

```text
notebooks/figures/
notebooks/supplementary_figures/
```

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


