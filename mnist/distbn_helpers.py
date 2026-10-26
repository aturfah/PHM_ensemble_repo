import scipy.stats as scistat
import numpy as np

from scipy.spatial import cKDTree
from scipy.special import gamma

def mix_sampler_1D(n, probs, distbns):
    output = []
    probs = probs / np.sum(probs)
    subsamples = scistat.multinomial(n, p=probs).rvs(1)[0]
    for idx in range(len(probs)):
        output.extend(distbns[idx].rvs(subsamples[idx]))

    return output

def mix_density_1D(x, probs, distbns):
    output = np.zeros(len(x))
    probs = probs / np.sum(probs)
    for idx in range(len(probs)):
        output += probs[idx] * distbns[idx].pdf(x)
    
    return output

def gmm_1d_generator(probs, locs, scales):
    assert len(probs) == len(locs)
    assert len(probs) == len(scales)
    M = len(probs)
    probs = probs / np.sum(probs)
    distbns = [scistat.norm(loc=locs[idx], scale=scales[idx]) for idx in range(M)]
    def density(x):
        return mix_density_1D(x, probs, distbns)

    def sampler(n):
        return mix_sampler_1D(n, probs, distbns)

    return density, sampler

def mvgmm_generator(probs, means, vars):
    assert len(probs) == len(means)
    assert len(probs) == len(vars)
    probs = probs / np.sum(probs)
    M = len(probs)
    def density(x):
        output = np.zeros(len(x))
        for idx in range(M):
            distbn = scistat.multivariate_normal(mean=means[idx], cov=vars[idx])
            output += probs[idx] * distbn.pdf(x)
        
        return output

    def sampler(n):
        D = len(means[0])
        output = np.zeros(shape=[0, D])
        subsamples = scistat.multinomial(n, p=probs).rvs(1)[0]
        for idx in range(M):
            distbn = scistat.multivariate_normal(mean=means[idx], cov=vars[idx])
            temp = distbn.rvs(subsamples[idx])
            if temp.shape[0] == 0:
                continue

            if subsamples[idx] == 1:
                temp = np.reshape(temp, (1, len(temp)))

            output = np.concatenate([output, temp], axis=0)

        return(output)

    return density, sampler


def mixt_1d_generator(probs, locs, dfs, scales):
    assert len(probs) == len(locs)
    assert len(probs) == len(dfs)
    assert len(probs) == len(scales)
    M = len(probs)
    distbns = [scistat.t(loc=locs[idx], df=dfs[idx], scale=scales[idx]) for idx in range(M)]

    def density(x):
        return mix_density_1D(x, probs, distbns)

    def sampler(n):
        return mix_sampler_1D(n, probs, distbns)

    return density, sampler


def knn_density(x, ref, k=5):
    """
    Computes the k-NN density estimation at points in 'x' based on 'ref'.
    
    Parameters:
    x (array-like): Input data where density is to be estimated, shape (n_queries, n_features).
    ref (array-like): Reference points forming the underlying distribution, shape (n_samples, n_features).
    k (int): The number of nearest neighbors to use.
    
    Returns:
    numpy.ndarray: Estimated density at each point in 'x'.
    """
    x = np.atleast_2d(x)
    ref = np.atleast_2d(ref)
    
    N, d = ref.shape
    
    # Build a KD-Tree from the reference points for fast nearest-neighbor lookups
    tree = cKDTree(ref)
    
    # Query the distance to the exact k-th nearest neighbor
    # tree.query with k=[k] returns only the k-th neighbor distance
    distances, _ = tree.query(x, k=[k])
    distances = distances.flatten()
    
    # Calculate the volume constant of a d-dimensional hypersphere
    # c_d = pi^(d/2) / Gamma(d/2 + 1)
    c_d = (np.pi ** (d / 2.0)) / gamma(d / 2.0 + 1.0)
    
    # Volume of the hypersphere defined by the radius (distance to k-th neighbor)
    volumes = c_d * (distances ** d)
    
    # Prevent division by zero in case of overlapping points
    volumes = np.maximum(volumes, 1e-15)
    
    # k-NN Density Estimate formula
    density = k / (N * volumes)
    
    return density