## Set directory for R installation
import os
os.environ["R_HOME"] = "/home/aturfah/.conda/envs/mnist_phm/lib/R"

## Disable BLAS/LAPACK multithreading to prevent multiprocessing deadlocks
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

from pathlib import path

## Imports
from functools import partial
from copy import deepcopy
from pathlib import Path
import plotnine as p9

import numpy as np
import pandas as pd

from sklearn.neighbors import KDTree
from mclustpy import mclustpy

## PHM Imports
from PHM import PHM
from PHMCluster import PHMCluster, deltaPmc
from distbn_helpers import mvgmm_generator

## Rpy2 imports for running R code
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

## File directories
file_directory = Path("/nfs/turbo/sph-xwen/aturfah/consensus/MNIST_test/data/")
raw_path = file_directory.joinpath("MNIST_raw_70k.npy")
labels_path = file_directory.joinpath("MNIST_labels_70k.npy")
pca_path = file_directory.joinpath("MNIST_pca_70k.npy")
tsne_path = file_directory.joinpath("MNIST_tsne_70k.npy")

save_dir = file_directory.joinpath("70k_subagg_gmm_tsne/")
part1_directory = file_directory.joinpath("70k_ref_models_tsne/")
part1_base = str(part1_directory.joinpath("{mthd}_tsne_repl{rpl}.json"))
save_dir = file_directory.joinpath("70k_subagg_gmm_tsne/")
save_dir_phm = file_directory.joinpath("70k_phm_replicates_tsne/")
phm_gmm_dir = file_directory.joinpath("70k_phm_gmm_replicates_tsne")
phm_skgmm_dir = file_directory.joinpath("70k_phm_skgmm_replicates_tsne")
phm_leiden_dir = file_directory.joinpath("70k_phm_leiden_replicates_tsne")
subagg_leiden_dir = file_directory.joinpath("70k_phm_subaggLeiden_replicates_tsne")
subagg_kmeans_dir = file_directory.joinpath("70k_phm_subaggKmeans_replicates_tsne")
phm_kmeans_dir = file_directory.joinpath("70k_phm_Kmeans_replicates_tsne")

if not part1_directory.exists():
    part1_directory.mkdir()

if not phm_leiden_dir.exists():
    phm_leiden_dir.mkdir()

if not phm_gmm_dir.exists():
    phm_gmm_dir.mkdir()

if not phm_skgmm_dir.exists():
    phm_skgmm_dir.mkdir()

if not save_dir.exists():
    save_dir.mkdir()

if not save_dir_phm.exists():
    save_dir_phm.mkdir()

if not subagg_leiden_dir.exists():
    subagg_leiden_dir.mkdir()

if not phm_kmeans_dir.exists():
    phm_kmeans_dir.mkdir()

if not subagg_kmeans_dir.exists():
    subagg_kmeans_dir.mkdir()

N = 70000
SUBSAMP_SIZE = 7000
REPLICATIONS = 100
CORES = 24
hdbscan_mcs = 80
leiden_resolution = 0.003575

## R imports
r_mclust = importr('mclust')
r_phm = importr('PHM')