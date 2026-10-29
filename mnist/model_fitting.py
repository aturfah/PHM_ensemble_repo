from pathlib import Path
import numpy as np
import sklearn.cluster as skcl
from sklearn.mixture import GaussianMixture as GMM
from sklearn.neighbors import KDTree
import hdbscan
import scanpy as sc
import sc3s

from PHMCluster import PHMCluster
from distbn_helpers import mvgmm_generator
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter


import dill
from pathos.multiprocessing import ProcessPool

from config import save_dir_mclust, r_phm

################################################
## Fit/Predict functions for baseline methods ##
################################################

def fit_predict_kmeans(X, n_clusters=10, random_state=0):
    return skcl.KMeans(n_clusters=n_clusters, random_state=random_state).fit_predict(X).tolist()


def fit_predict_hclust(X, n_clusters=10, random_state=0):
    return skcl.AgglomerativeClustering(n_clusters=n_clusters).fit_predict(X).tolist()


def fit_predict_skgmm(X, n_clusters=10, random_state=0, init="random"):
    return GMM(n_components=n_clusters, 
               random_state=random_state,
               max_iter=500, 
               init_params=init).fit_predict(X)


def fit_predict_mclust2(X, n_clusters=10, random_state=0):
    """
    set.seed(idx)
    mcl <- Mclust(X, G=10)
    """
    ## Load the R file
    fname = save_dir_mclust.joinpath("mclustG{}_pred_r{}.npy".format(n_clusters, random_state))
    if not fname.exists():
        raise FileNotFoundError(fname)

    return np.load(fname)


def fit_predict_sc3s(X, n_clusters=10, random_state=0, n_runs=5):
    adata = sc.AnnData(X)
    adata.obsm["X_pca"] = X
    sc3s.tl.consensus(adata, 
                      n_clusters=n_clusters, 
                      d_range=X.shape[1], 
                      n_runs=n_runs, 
                      random_state = random_state)
    return adata.obs['sc3s_{}'.format(n_clusters)].tolist()


def fit_predict_leiden(X, resolution=0.0007, random_state=0):
    adata = sc.AnnData(X)
    adata.obsm["X_pca"] = X
    sc.pp.neighbors(adata, random_state=random_state)
    sc.tl.leiden(adata, random_state=random_state, resolution=resolution)
    return adata.obs['leiden'].tolist()


def fit_predict_hdbscan(X, mcs=30, **kwargs):
    hdbscan_clust = hdbscan.HDBSCAN(min_cluster_size=mcs).fit_predict(X)
    adjusted_preds = np.copy(hdbscan_clust)
    ## Turn noise points into singletons
    noise_mask = adjusted_preds == -1
    num_noise = np.sum(noise_mask)
    adjusted_preds[noise_mask] = X.shape[0] + np.arange(num_noise)
    return adjusted_preds.tolist()

#######################################
## Convert R params object to Python ##
#######################################
def convert_R_params_object2(res, probscale=1):
    clusters = []
    densities = []
    with localconverter(robjects.default_converter + numpy2ri.converter):
        for clust_params in res:
            ## Extract params into dictionary numpy arrays
            param_dict = {}
            for val in clust_params.items():
                param_dict[val.name] = val.value

            ## Construct PHMClusters for each mixture component
            G = param_dict["prob"].shape[0]
            for comp_idx in range(G):
                prob = param_dict["prob"][comp_idx]
                mean = np.reshape(param_dict["mean"][:, comp_idx], (param_dict["mean"].shape[0], ))
                var = np.reshape(param_dict["var"][:, :, comp_idx], (param_dict["var"].shape[0], param_dict["var"].shape[1]))
                density, sampler = mvgmm_generator([prob * probscale], [mean], [var])
                clust = PHMCluster([prob * probscale], [sampler], [density])
                clusters.append(clust)

    return clusters, densities

##########################################
## PHM Density Estimation and Retrieval ##
##########################################
def subaggDensityEstimation(X, replicates, subSampleSize, fit_func, extract_func,
                            replicate_seeds = None,
                            numCores=1,
                            save_directory=None, prefix=None,
                            **kwargs):
    save_path = Path(save_directory)
    if not save_path.exists():
        save_path.mkdir()
    if replicate_seeds is None:
        replicate_seeds = [x for x in range(replicates)]
    
    def _replicate_func(repl, seed=None):
        if seed is None:
            seed = repl

        save_fname = None
        
        ## If saving, try to load
        if save_directory is not None:
            fname = "{}_repl{}.pkl".format(prefix, repl)
            save_fname = save_path.joinpath(fname)
            if save_fname.exists():
                ## load dill
                with save_fname.open("rb") as sf:
                    dens = dill.load(sf)
                    total_prob = np.sum([x.alpha() for x in dens])
                    for x in dens:
                        x.probs /= total_prob * replicates
                    return dens

        np.random.seed(seed)
        N = X.shape[0]
        sub_X = X[np.random.randint(low=0, high=N, size=subSampleSize), ]
        res_raw = fit_func(sub_X, rs=seed, **kwargs)
        res = extract_func(res_raw, 1/replicates)

        ## Write if we are in fact saving
        if save_fname is not None:
            ## write dill
            with save_fname.open("wb") as sf:
                dill.dump(res, sf)

        return res

    output = []
    
    if numCores != 1:
        pool = ProcessPool(numCores)
        with pool:
            for res in pool.uimap(_replicate_func, range(replicates), replicate_seeds):
                output.extend(res)
    else:
        for repl in range(replicates):
            res = _replicate_func(repl)
            output.extend(res)

    return output

## Fit/Extract Functions for Leiden and GMM
def skGMM_ff(X, G_range, rs, **kwargs):
    ## Fit candidate models
    np.random.seed(rs)
    temp = []
    for G in G_range:
        temp.append(GMM(n_components=G, **kwargs).fit(X))

    ## Model with lowest BIC
    bic_vec = [x.bic(X) for x in temp]
    idx = np.argmin(bic_vec)
    return temp[idx]

def skGMM_extract(res, wgt_scale):
    clusters = []
    densities = []
    for idx in range(res.n_components):
        cov = res.covariances_[idx]
        mean = res.means_[idx]
        wgt = res.weights_[idx]
    
        density, sampler = mvgmm_generator([wgt], [mean], [cov])
        densities.append(density)
        clusters.append(
            PHMCluster([wgt * wgt_scale], [sampler], [density])
        )
    
    return clusters

def leiden_ff(X, rs=1, resolution=0.05, weighting_neighbors=1, **kwargs):
    ## Run Leiden clustering
    np.random.seed(rs)
    pred_leiden = fit_predict_leiden(X, resolution=resolution)
    pred_leiden = np.asarray([int(x) for x in pred_leiden])

    cluster_order = []
    for x in pred_leiden:
        if x in cluster_order:
            continue

        cluster_order.append(x)
    
    ## Compute weights for weighted density estimation
    cluster_distance = np.zeros(shape=(X.shape[0], len(cluster_order)))
    for label_idx in range(len(cluster_order)):
        ## Get observations in cluster
        label = cluster_order[label_idx]
        label_mask = np.array([True if x == label else False for x in pred_leiden]) 
        X_label = X[label_mask, ]
    
        ## Tree for nearest neighbor
        tree = KDTree(X_label)
        k_query = 2 if weighting_neighbors == 1 else weighting_neighbors 
        dist, ind = tree.query(X, k=k_query)
    
        if weighting_neighbors == 1:
            ## Nearest non-self neighbor in cluster
            dist = np.where(label_mask, dist[:, 1], dist[:, 0])
        else:
            dist = dist[:, weighting_neighbors - 1]
            
        cluster_distance[:, label_idx] = dist

    cluster_weight = np.exp(-cluster_distance)
    cluster_weight = cluster_weight / cluster_weight.sum(axis=1, keepdims=True)

    ## Weighted Density Estimate
    with localconverter(robjects.default_converter + numpy2ri.converter):
        weighted_partition_density = r_phm.constructPmcParamsWeightedPartition(
            np.asarray(pred_leiden),
            X,
            weights=cluster_weight,
            **kwargs
        )

    return weighted_partition_density

def leiden_extract(res, wgt_scale):
    ## Convert to form to be consumed
    split_clusters, _ = convert_R_params_object2(res, wgt_scale)

    return split_clusters