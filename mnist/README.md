# MNIST Analysis

Code provided to replicate the MNIST analysis. Generally to run please execute scripts in the following order.

**Make sure to set the `R_HOME` environment variable in `config.py`**

```
## Install necessary packages
Rscript mnist_density_estimation.R
python run_cluster_performance.py
python run_phm_visualizations.py 
```

## Environment Set-Up

Below are the commands to set up the conda environment to run this analysis.

```
conda create --name PHM_ensemble python=3.11 -y

## If running with Jupyter notebook, set an environment for jupyter to use
conda install -n PHM_ensemble ipykernel -y
conda run -n PHM_ensemble python -m ipykernel install --user --name=PHM_ensemble --display-name "Python (PHM Ensemble)"

## Install R and necessary packages
conda install -n PHM_ensemble conda-forge::r conda-forge::r-base -y
conda install -n PHM_ensemble -c conda-forge r-tidyverse r-mclust r-aricode r-reticulate r-clue r-fastcluster r-remotes r-biocmanager

conda activate PHM_ensemble
Rscript install_packages.R

## Install necessary Python Packages
conda run -n PHM_ensemble pip install pathos dill plotnine pandas scikit-learn opentsne scanpy rpy2 torchvision sc3s mclustpy igraph leidenalg hdbscan
```

## File Information

Below is a description of the files provided grouped by functionality.


### PHM Python implementation

For ease of use this is an initial Python implementation of the [`PHM` R package](https://github.com/aturfah/PHM). 

- `distbn_helpers.py`: Implements basic functions to generate Gaussian mixture PDFs and samplers using `scipy`.
- `PHMCluster.py`: Implements a general cluster mixture density class (functions analogous to the scipy library, ex: `pdf()` and `rvs()`) with capability to merge with another cluster yielding a combined density. Additionally implements the functions to evaluate Pmc and compute the delta Pmc matrix.
- `PHM.py`: Implements the PHM algorithm and associated visualizations. Also allows for evaluations along the merging process, e.g., compute the Posterior matrix at a given value of *K*.


### General analysis functions and settings

- `config.py`: Necessary imports/configuration settings, prepares file directories for different outputs.
- `prepare_data.py`: Implements functions to download and preprocess the MNIST data.
- `helpers.py`: Implements json saving/loading functions, multiprocess-safe printing, as well as functions for initializing/calculating performance.

### Comparison to baseline methods

- `run_cluster_performance.py`: Implements the loops to evaluate and compare the PHM variants (single-replicate, localized ensemble, global ensemble) with the baseline methods. Outputs a table via text and a boxplot for cluster performance.
- `mnist_density_estimation.R`: Performs density estimation for, PHM (single replicate), as well as PHM-Ensemble (global). Also runs clustering for mclust, hierarchical clustering, as well as the Oracle model.
- `baseline.py`: Implements function to evaluate the baseline methods, stores clustering results.
- `model_fitting.py`: Implements basic functions to cluster using the baseline methods (`fit_predict_*`), as well as the localized density estimation using Leiden clustering.
- `phm_single.py`, `phm_ensemble_global.py`, `phm_ensemble_localized.py`: Implements the functions to run the PHM algorithm on the corresponding density estimate. Generates cluster output.


### PHM dendrogram and subclustering  analysis

- `run_phm_visualizations`: Loads the first PHM-ensemble (global) replicate and generates the corresponding dendrogram/heatmap visualization. Also generates drill-down image for subcluster analysis.
- `vis_helpers.py`: Implements the digit-grid visualization for subcluster analysis.
