# Human Immune Health Atlas Analysis Scripts

Code provided to replicate the Human Immune Health Atlas (HIHA) analysis. Generally to run please execute scripts in the following order. Please see the MNIST folder for information on setting up the environment.

```
## Install necessary packages
sh setup_symlinks.sh
Rscript density_estimation.R
python profile_phm.py
python generate_visualizations.py
```

Please download the full HIHA scRNA dataset from the [Allen Institute Website](https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/scrna/).

## File Descriptions

- `setup_symlinks.sh`: Creates symlinks to the `PHM.py`, `PHMCluster.py`, and `distbn_helpers.py` files from the MNIST folder for use here.
- `density_estimation.R`: Runs global ensemble density esstimation on the HIHA dataset
- `profile_phm.py`: Runs the PHM algorithm based on the ensemble density estimate and produces time/memory usage estimates.
- `generate_visualizations.py`: Produces cell-type and T/NK cell drill-down visualizations from the paper.
- `config.py`: Points to R environment folders and sets up filenames.
- `helpers.py`: Functions for loading R objects and constructing contingency tables based on the posterior matrix for the dendrogram visualizations.