## Set directory for R installation
import os
## Example: /home/aturfah/.conda/envs/PHM_ensemble/lib/R/
os.environ["R_HOME"] = "R/environment/location"

## Disable BLAS/LAPACK multithreading to prevent multiprocessing deadlocks
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

## Imports
from pathlib import Path
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

## R imports
r_mclust = importr('mclust')
r_phm = importr('PHM')

## File directories
file_directory = Path("output/")
images_directory = file_directory.joinpath("images")
raw_path = file_directory.joinpath("MNIST_raw_70k.npy")
labels_path = file_directory.joinpath("MNIST_labels_70k.npy")
pca_path = file_directory.joinpath("MNIST_pca_70k.npy")
tsne_path = file_directory.joinpath("MNIST_tsne_70k.npy")

baseline_directory = file_directory.joinpath("70k_baseline/")
part1_base = str(baseline_directory.joinpath("{mthd}_tsne_repl{rpl}.json"))
save_dir_global_ensemble = file_directory.joinpath("70k_ensemble_global/")
save_dir_phm_single = file_directory.joinpath("70k_single_PHM/")
save_dir_mclust = file_directory.joinpath("70k_mclustG10")
save_dir_localized_ensemble = file_directory.joinpath("70k_ensemble_localized/")

if not file_directory.exists():
    file_directory.mkdir()

if not images_directory.exists():
    images_directory.mkdir()

if not baseline_directory.exists():
    baseline_directory.mkdir()

if not save_dir_phm_single.exists():
    save_dir_phm_single.mkdir()

if not save_dir_global_ensemble.exists():
    save_dir_global_ensemble.mkdir()

if not save_dir_localized_ensemble.exists():
    save_dir_localized_ensemble.mkdir()

if not save_dir_mclust.exists():
    save_dir_mclust.mkdir()

# if not phm_leiden_dir.exists():
#     phm_leiden_dir.mkdir()

# if not phm_gmm_dir.exists():
#     phm_gmm_dir.mkdir()

# if not phm_skgmm_dir.exists():
#     phm_skgmm_dir.mkdir()

# if not save_dir.exists():
#     save_dir.mkdir()

# if not subagg_leiden_dir.exists():
#     subagg_leiden_dir.mkdir()

# if not phm_kmeans_dir.exists():
#     phm_kmeans_dir.mkdir()

# if not subagg_kmeans_dir.exists():
#     subagg_kmeans_dir.mkdir()


## Parameters/other things for rest of files
N = 70000
SUBSAMP_SIZE = 7000
# REPLICATIONS = 10
REPLICATIONS = 100
CORES = 30
hdbscan_mcs = 80
leiden_resolution = 0.003575