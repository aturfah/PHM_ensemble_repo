import numpy as np
import scipy.stats as scistat
from typing import List, Callable


class PHMCluster:
    def __init__(self, alpha:np.array, sampler:List[Callable], density:List[Callable], sampler_distbn:List[Callable]=None):
        if not isinstance(alpha, np.ndarray):
            alpha = np.array(alpha)

        self.K = alpha.shape[0]
        self.probs = alpha
        self.pdf_funcs = density
        self.rvs_funcs = sampler
        self.spdf_funcs = sampler_distbn

    def rvs(self, n:int):
        sub_probs = self.probs / np.sum(self.probs)
        group_counts = scistat.multinomial(n, sub_probs).rvs(1)[0]

        output = None
        for idx in range(len(group_counts)):
            gc = group_counts[idx]
            if gc == 0:
                continue

            samp = self.rvs_funcs[idx](gc)
            if output is None:
                output = samp
            else:
                ## Zero-dimensional array
                output = np.concatenate([output, samp], axis=0)

        return output

    def pdf(self, x:np.ndarray, func:str="pdf"):
        sub_probs = self.probs / np.sum(self.probs)
        if func == "pdf":
            func_list = self.pdf_funcs
        elif func == "spdf":
            func_list = self.spdf_funcs
        else:
            raise ValueError("func must be either 'pdf' or 'spdf'")

        pdf_mat = np.concatenate([func(x) for func in func_list]).reshape((-1, len(func_list)), order="F")
        temp = np.multiply(pdf_mat, sub_probs)

        return np.sum(temp, axis=1)

    def log_sampling_ratio(self, x:np.ndarray):
        if self.spdf_funcs is None:
            ## Assume that density and sampler density are the same
            out = np.zeros(x.shape[0]) - np.log(x.shape[0])
            return (out, out)

        num = np.log(self.pdf(x, "pdf"))
        den = np.log(self.pdf(x, "spdf"))
        return (num, den)


    def alpha(self):
        return(np.sum(self.probs))

    def merge(self, obj):
        self.probs = np.concatenate([self.probs, obj.probs], 0)
        self.rvs_funcs.extend(obj.rvs_funcs)
        self.pdf_funcs.extend(obj.pdf_funcs)

        ## Sampler pdf functions
        if self.spdf_funcs is None and obj.spdf_funcs is None:
            return
        else:
            if self.spdf_funcs is None:
                self.spdf_funcs = [lambda x: np.zeros(x.shape[0]) + 1/self.K for _ in range(self.K)]
            
            if obj.spdf_funcs is None:
                spdf_obj = [lambda x: np.zeros(x.shape[0]) + 1/obj.K for _ in range(obj.K)]
            else:
                spdf_obj = obj.spdf_funcs

            self.spdf_funcs.extend(spdf_obj)

def sampleMixture(n:int, distbn:List[PHMCluster]) -> np.ndarray:
    output = None
    weight_vec = [x.alpha() for x in distbn]

    comp_sizes = scistat.multinomial(n, weight_vec).rvs(1)[0]
    for idx in range(len(comp_sizes)):
        csize = comp_sizes[idx]
        new_samp = distbn[idx].rvs(csize)
        if new_samp is None:
            continue
        elif output is None:
            output = new_samp
        else:
            output = np.concatenate([output, new_samp])

    return output

def computePosterior(mat:np.ndarray, distbn:List[PHMCluster]) -> np.ndarray:
    K = len(distbn)
    dens_mat = np.zeros(shape=[mat.shape[0], K])
    for idx in range(K):
        dens_mat[:, idx] = distbn[idx].alpha() * distbn[idx].pdf(mat)

    norm_mat = np.sum(dens_mat, axis=1)
    return dens_mat / norm_mat[:, None]


def deltaPmc(distbn:List[PHMCluster], 
                    mcSamples:int=1e5, 
                    batchSize:int=1e4,
                    seed=1) -> np.ndarray:
    np.random.seed(seed)

    sample_sizes = [batchSize for _ in range(int(np.ceil(mcSamples/batchSize)))]
    if np.sum(sample_sizes) > mcSamples:
        sample_sizes[-1] -= int(np.sum(sample_sizes) - mcSamples)

    K = len(distbn)
    output = np.zeros((K, K))
    for ssize in sample_sizes:
        sample = sampleMixture(ssize, distbn)

        post_mat = computePosterior(sample, distbn)
        temp = (post_mat.T @ post_mat) / mcSamples
        np.fill_diagonal(temp, 0)

        output += temp
        # for idx1 in range(K):
        #     for idx2 in range(idx1+1, K):
        #         val = np.sum(
        #             post_mat[:, idx1] * post_mat[:, idx2]
        #         ) / mcSamples
        #         output[idx1, idx2] += val
        #         output[idx2, idx1] += val

    return output

def pmc(distbn:List[PHMCluster], 
                mcSamples:int=1e5, 
                batchSize:int=1e4,
                seed=1) -> float:
    np.random.seed(seed)

    sample_sizes = [batchSize for _ in range(int(np.ceil(mcSamples/batchSize)))]
    if np.sum(sample_sizes) > mcSamples:
        sample_sizes[-1] -= int(np.sum(sample_sizes) - mcSamples)

    K = len(distbn)
    output = 0
    for ssize in sample_sizes:
        sample = sampleMixture(ssize, distbn)

        post_mat = computePosterior(sample, distbn)
        output += np.sum(post_mat * (1 - post_mat)) / mcSamples

    return output

if __name__ == "__main__":
    distbn1 = scistat.norm(loc=0, scale=1)
    distbn2 = scistat.norm(loc=2, scale=1)
    pdf1 = lambda x: distbn1.pdf(x)
    pdf2 = lambda x: distbn2.pdf(x)
    rvs1 = lambda x: distbn1.rvs(x)
    rvs2 = lambda x: distbn2.rvs(x)

    comp = PHMCluster([0.2, 0.4], sampler=[rvs1, rvs2], density=[pdf1, pdf2])

    samp = comp.rvs(50)
    dens = comp.pdf(samp)

    samp = comp.rvs(1)
    dens = comp.pdf(samp)

    comp2 = PHMCluster([0.4], sampler=[rvs1], density=[pdf1])
    # comp.merge(comp2)

    ## Pmc with this class
    print(deltaPmc([comp, comp2]))
    print(pmc([comp, comp2]))