import numpy as np
from PHMCluster import deltaPmc, PHMCluster, pmc, computePosterior
from typing import List, Callable
from copy import deepcopy
import scipy
import pandas as pd
from itertools import product

import plotnine as p9

from pprint import pprint

class PHM:
    def __init__(self, distbn:List[PHMCluster], scaling):
        self.mat = None
        ## Distribution information
        self.c_prob = [x.alpha() for x in distbn]
        self.K = len(self.c_prob)
        self.distbn = distbn
        self.scalePHM = scaling

        ## PHM Tracking
        self.results = False
        self.mergeComponents = []
        self.mergeNames = []
        self.mergeValues = []
        self.mergeDeltaPmc = []

        ## Die
        self.plot_df = None
        self.order_x = None


    def computeDeltaPmc(self, override=False, **kwargs):
        if self.mat is not None and not override:
            ## Delta Pmc has already been calculated
            raise RuntimeError("Delta Pmc matrix is already computed. Set override=True to recompute.")

        self.mat = deltaPmc(self.distbn, **kwargs)

        return self.mat

    def run(self, override=False, **kwargs):
        ## Compute Delta Pmc if not already there
        if self.mat is None:
            self.computeDeltaPmc(**kwargs)

        if self.results and not override:
            raise RuntimeError("PHM has already ben run. set override=False to recompute.")

        ## Clear previous results
        self.mergeComponents = []
        self.mergeNames = []
        self.mergeValues = []
        self.mergeDeltaPmc = []

        ## Tracking for clusters
        cluster_sizes = deepcopy(self.c_prob)
        num_components = [1 for _ in range(self.K)]
        cluster_names = [x for x in range(self.K)]

        ## Prepare matrix for PHM
        temp = deepcopy(self.mat)
        np.fill_diagonal(temp, -1)

        for _ in range(self.K-1):
            if self.scalePHM == "alpha":
                value_matrix = np.multiply(temp, 1 / np.outer(cluster_sizes, cluster_sizes))
            elif self.scalePHM == "average":
                value_matrix = np.multiply(temp, 1 / np.outer(num_components, num_components))
            elif self.scalePHM == "unscaled":
                value_matrix = temp
            else:
                raise ValueError("Unsupported scaling: {}".format(self.scalePHM))

            merge_idx = np.unravel_index(np.argmax(value_matrix),
                                         value_matrix.shape)
            ## Ensure we're going in upper triangle
            if merge_idx[0] > merge_idx[1]:
                merge_idx = [merge_idx[1], merge_idx[0]]

            ## Track this
            max_val = value_matrix[merge_idx]
            self.mergeComponents.append(merge_idx)
            self.mergeValues.append(max_val)
            self.mergeDeltaPmc.append(temp[merge_idx])
            self.mergeNames.append([cluster_names[merge_idx[0]], 
                                    cluster_names[merge_idx[1]]])

            ## Update matrix
            temp[:, merge_idx[0]] += temp[:, merge_idx[1]]
            temp[merge_idx[0], :] += temp[merge_idx[1], :]
            temp[merge_idx[0], merge_idx[0]] = -1

            ## Zero out elements
            temp[:, merge_idx[1]] = 0
            temp[merge_idx[1], :] = 0

            ## Update cluster sizes/counts
            num_components[merge_idx[0]] += num_components[merge_idx[1]]
            cluster_sizes[merge_idx[0]] += cluster_sizes[merge_idx[1]]
            cluster_names[merge_idx[0]] = "{}|{}".format(cluster_names[merge_idx[0]],
                                                         cluster_names[merge_idx[1]])

        self.results = True
        self._postprocess()

    def runv2(self, override=False, **kwargs):
        ## Compute Delta Pmc if not already there
        if self.mat is None:
            self.computeDeltaPmc(**kwargs)

        if self.results and not override:
            raise RuntimeError("PHM has already ben run. set override=False to recompute.")

        ## Clear previous results
        self.mergeComponents = []
        self.mergeNames = []
        self.mergeValues = []
        self.mergeDeltaPmc = []

        ## Tracking for clusters
        cluster_sizes = deepcopy(self.c_prob)
        num_components = [1 for _ in range(self.K)]
        cluster_names = [x for x in range(self.K)]

        ## Prepare matrix for PHM; lower triangle is all -1
        temp = 2 * deepcopy(self.mat)
        np.fill_diagonal(temp, -1)

        while temp.shape != (1, 1):
            if self.scalePHM == "alpha":
                value_matrix = np.multiply(temp, 1 / np.outer(cluster_sizes, cluster_sizes))
            elif self.scalePHM == "average":
                value_matrix = np.multiply(temp, 1 / np.outer(num_components, num_components))
            elif self.scalePHM == "unscaled":
                value_matrix = temp
            else:
                raise ValueError("Unsupported scaling: {}".format(self.scalePHM))

            merge_idx = np.unravel_index(np.argmax(value_matrix),
                                                    value_matrix.shape)
            ## Ensure we're going in upper triangle
            if merge_idx[0] > merge_idx[1]:
                merge_idx = [merge_idx[1], merge_idx[0]]

            if merge_idx[0] >= merge_idx[1]:
                raise RuntimeError("Maximum value found in lower triangle or along diagonal.")

            ## Track this
            max_val = value_matrix[merge_idx]
            self.mergeComponents.append(merge_idx)
            self.mergeValues.append(max_val)
            self.mergeDeltaPmc.append(temp[merge_idx])
            self.mergeNames.append([cluster_names[merge_idx[0]], 
                                    cluster_names[merge_idx[1]]])

            ## Update matrix
            temp[:, merge_idx[0]] += temp[:, merge_idx[1]]
            temp[merge_idx[0], :] += temp[merge_idx[1], :]
            temp[merge_idx[0], merge_idx[0]] = -1

            ## Update cluster sizes/counts
            num_components[merge_idx[0]] += num_components[merge_idx[1]]
            cluster_sizes[merge_idx[0]] += cluster_sizes[merge_idx[1]]
            cluster_names[merge_idx[0]] = "{}|{}".format(cluster_names[merge_idx[0]],
                                                         cluster_names[merge_idx[1]])

            ## Remove unnecesssary elements
            temp = np.delete(temp, merge_idx[1], 0)
            temp = np.delete(temp, merge_idx[1], 1)
            del num_components[merge_idx[1]]
            del cluster_sizes[merge_idx[1]]
            del cluster_names[merge_idx[1]]

        self.results = True
        self._postprocess()

    def _postprocess(self):
        if not self.results:
            raise RuntimeError("Please run PHM.run() before calling this (though you shouldn't be calling this directly)")

        ## Originally in order of merges (idx0 is first merge, idkK is K+1th merge)
        ## Change this so that idxK has K+1 clusters left because that is how code was originally written
        self.mergeComponents = list(reversed(self.mergeComponents))
        self.mergeNames = list(reversed(self.mergeNames))
        self.mergeValues = list(reversed(self.mergeValues))
        self.mergeDeltaPmc = list(reversed(self.mergeDeltaPmc))

    def _prepareVisData(self, initK=None, scale_heights="log10", group_probs=None, threshold=0):
        """Prepare data for visualization"""
        if initK is None:
            initK = self.K

        if initK > self.K:
            raise ValueError("Chosen value of initK larger than original number of clusters")

        ## Track whether we want to visualize group probabilities across PHM dend
        gprobs_unspec = False
        if group_probs is None:
            gprobs_unspec = True
            group_probs = [1 for _ in range(initK)]

        below_threshold = [x <= threshold for x in self.mergeValues]
        no_valid_merges = np.all(below_threshold)
        all_valid_merges = np.all([not x for x in below_threshold])
        stop_idx = None
        if all_valid_merges:
            stop_idx = -1
        elif not no_valid_merges:
            stop_idx = np.argmax([x <= threshold for x in self.mergeValues])

        ## Determine the height to be used
        if no_valid_merges:
            ## Everything is at same height; 1
            height = 1 + np.zeros((initK))
        else:
            height = self.mergeValues
            if self.scalePHM == "unscaled":
                clust_sizes = [1 for _ in range(self.K)]
                sizes = np.zeros((self.K-1))

                for idx in range(self.K-1, 0, -1):
                    idx = idx - 1
                    mc = self.mergeComponents[idx]
                    mcL = mc[0]
                    mcR = mc[1]
                    clust_sizes[mcL] += clust_sizes[mcR]
                    del clust_sizes[mcR]
                    sizes[idx] = clust_sizes[mcL]

                ## Scale by number of combinations
                height = np.array(height)
                height = height / scipy.special.comb(sizes, 2)

            ## Only get the valid heights
            height = height[(stop_idx+1):initK]


        ## Scale the height based on the merging rule
        if scale_heights == "unscaled":
            height = 1 / height
        elif scale_heights == "log10":
            if self.scalePHM == "alpha":
                ## Scale so that height is interpretable?
                height /= np.max(height)

            height = -1 * np.log(height)
        elif scale_heights == "log10":
            height = -1 * np.log10(height)

        ## Slight height offset so that things are guaranteed to not be at same height
        height += np.array([x for x in range(len(height))]) * 1e-6

        ## Construct dataframe with height information for merges
        output = []
        group_probs_new = deepcopy(group_probs)
        alpha_vec = self.recoverParams(initK, ["prob"])["prob"]
        merge_tree = [x for x in range(initK)]
        component_id_map = [x for x in range(initK)]
        height_tracker = [{"base": 0, "height": 0} for _ in range(initK)]
        base_height = 0
        if no_valid_merges:
            for idx in range(initK):
                output.append({
                    "ID": idx,
                    "y": 0,
                    "yend": 1 + idx * 1e-6,
                    "gprob": group_probs_new[idx],
                    "linetype": "dashed"
                })
        else:
            merge_components = self.mergeComponents[(stop_idx+1):(initK-1)]
            for idx in reversed(range(len(merge_components))):
                mc = merge_components[idx]
                hgt = height[idx]

                ## Make sure height is non-decreasing along this tree
                for idx2 in mc:
                    hgt = np.max([hgt, height_tracker[idx2]["base"]])

                for idx2 in mc:
                    height_tracker[idx2]["height"] = hgt

                for idx2 in mc:
                    output.append({
                        "ID": component_id_map[idx2],
                        "y": height_tracker[idx2]["base"],
                        "yend": height_tracker[idx2]["height"],
                        "gprob": group_probs_new[idx2],
                        "linetype": "solid"
                    })

                ## TODO: Update group probs
                if gprobs_unspec:
                    group_probs_new = group_probs_new[1:len(group_probs)]
                else:
                    prob1 = alpha_vec[mc[0]]
                    prob2 = alpha_vec[mc[1]]
                    group_probs_new[mc[0]] = (group_probs_new[mc[0]] * prob1 \
                        + group_probs_new[mc[1]] * prob2) / \
                        (prob1 + prob2)
                    del group_probs_new[mc[1]]
                    alpha_vec[mc[0]] += alpha_vec[mc[1]]
                    np.delete(alpha_vec, mc[1])

                ## Update height tracker
                base_height = hgt
                height_tracker[mc[0]]["base"] = base_height
                del height_tracker[mc[1]]

                ## Update Component ID Map
                component_id_map.remove(component_id_map[mc[1]])

                ## Update Merge Tree
                merge_tree[mc[0]] = {
                    "left": merge_tree[mc[0]],
                    "right": merge_tree[mc[1]],
                    "order": idx
                }
                del merge_tree[mc[1]]


        ## Get ordering of observations along X axis
        def order_x_axis(merge_res):
            if isinstance(merge_res, dict):
                return [order_x_axis(merge_res["left"]),
                        order_x_axis(merge_res["right"])]
            else:
                return merge_res

        def flatten(nest_list):
            for nl in nest_list:
                if isinstance(nl, list):
                    yield from flatten(nl)
                else:
                    yield nl

        if len(merge_tree) == 1:
            order_x = order_x_axis(merge_tree[0])
        else:
            order_x = [order_x_axis(mt) for mt in merge_tree]

        order_x = list(flatten(order_x))

        ## Premature stop means we may have unmerged clusters
        if not all_valid_merges:
            max_height = np.max([x["yend"] for x in output])
            unmerged_ids = [x for x in merge_tree if not isinstance(x, dict)]
            output_ids = [x["ID"] for x in output]
            ## Accounts for bypass where no_valid_merges
            unmerged_ids = set(unmerged_ids).difference(output_ids)
            for uid in unmerged_ids:
                output.append({
                    "ID": uid,
                    "y": 0,
                    "yend": max_height + (1+uid) * (1e-6 * max_height),
                    "gprob": group_probs[uid],
                    "linetype": "dashed"
                })

        ## Assign X positions to base level merges
        map_xposns = lambda val: order_x.index(val)
        for obj in output:
            if obj["y"] == 0:
                obj["x"] = map_xposns(obj["ID"])
            else:
                obj["x"] = None

            obj["xend"] = None

        df = pd.DataFrame(output)

        ## Assign X positions to merges further up trees
        ## Subtree yend corresponds to these branches' y value; must be unique
        while np.any(np.isnan(df["x"])):
            merged = pd.merge(
                df, 
                df, 
                left_on='y', 
                right_on='yend', 
                how='left', 
                suffixes=('', '.y')
            )
            group_cols = ['ID', 'y', 'yend', 'x', 'gprob', 'linetype']
            agg_df = merged.groupby(group_cols, dropna=False)['x.y'].agg(lambda s: s.mean(skipna=False)).reset_index()
            agg_df['xend'] = np.where(agg_df['x'].isna(), agg_df['x.y'], agg_df['x'])
            agg_df['x'] = agg_df['xend']
            agg_df = agg_df.sort_values(by='yend', ascending=False)
            df = agg_df.drop(columns=['x.y']).reset_index(drop=True)

        ## Add Horizontal Components
        horiz_comps = df.groupby(['yend'], dropna=False).agg(
            xend=('x', 'max'),
            x=('x', 'min'),
            gprob=('gprob', 'mean')
        ).reset_index()
        horiz_comps['y'] = horiz_comps['yend']
        horiz_comps['linetype'] = 'solid'

        df = pd.concat([df, horiz_comps])

        ## Finally make sure vertical lines are accounted for
        df['xend'] = df['xend'].infer_objects(copy=False).fillna(df['x'])

        self.plot_df = df
        self.order_x = order_x

    def dendrogram(self,
                   initK:int=None,
                   scale_heights: str="log10",
                   group_probs:np.array=None,
                   threshold:float=0, 
                   groupColorMax:str="#000000",
                   groupColorMin:str="#CCCCCC") -> p9.ggplot:
        self._prepareVisData(initK, scale_heights=scale_heights, group_probs=group_probs, threshold=threshold)
        if initK is None:
            initK = self.K

        plt = p9.ggplot(self.plot_df, p9.aes(x="x", y="y", xend="xend", yend="yend")) + \
            p9.geom_segment(p9.aes(linetype="linetype", color="gprob")) + \
            p9.scale_linetype_identity() + \
            p9.scale_color_gradient(low=groupColorMin, high=groupColorMax) + \
            p9.xlab("") + p9.ylab("") + \
            p9.theme_bw() + \
            p9.scale_x_continuous(
                breaks=list(range(0, initK)), 
                limits=(0, initK),
                labels=self.order_x) + \
            p9.scale_y_continuous(expand=(0, 0, 0.05, 0)) + \
            p9.theme(
                # axis_text_x=displayAxisFmt, 
                axis_text_y=p9.element_blank(),
                axis_ticks_y=p9.element_blank(),
                panel_grid=p9.element_blank(),
                panel_spacing=0, 
                panel_border=p9.element_blank(),
                legend_position="none"
            )

        return plt

    def heatmap(self, 
                initK:int=None, 
                scale_heights:str="log10", 
                group_probs:np.array=None, 
                threshold:float=0, 
                plot_lims:tuple = None,
                phm:bool=True,
                scale_axis=False):
        self._prepareVisData(initK, scale_heights=scale_heights, group_probs=group_probs, threshold=threshold)
        if initK is None:
            initK = self.K


        if phm is True:
            mat = np.zeros((initK, initK))
            label_map = [
                [x] for x in range(initK)
            ]
            for idx in range(initK-2, -1, -1):
                mc = self.mergeComponents[idx]
                mv = self.mergeValues[idx]
                labels_0 = label_map[mc[0]]
                labels_1 = label_map[mc[1]]

                mat_indices = [x for x in product(labels_0, labels_1)]
                for ind in mat_indices:
                    mat[ind[0], ind[1]] = mv
                    mat[ind[1], ind[0]] = mv

                label_map[mc[0]].extend(label_map[mc[1]])
                del label_map[mc[1]]
        else:
            mat = self.recoverDeltaPmc(initK)
            # if self.scalePHM == "alpha":
            #     alpha_vec = self.recoverParams(initK)["prob"]
            #     mat = np.multiply(mat, np.outer(alpha_vec, alpha_vec))

        mat_long = (pd.DataFrame(mat)
                    .assign(X=lambda df: range(len(df)))
                    .melt(id_vars="X", var_name="Y", value_name="Z")
                    .assign(Z_old=lambda df: df["Z"],
                            Z=lambda df: np.log10(df["Z"].where(df["Z"] > 0)))
                    )

        if plot_lims is None:
            plot_lims = [np.min(mat_long["Z"]), np.max(mat_long["Z"])]

        mid_point = (plot_lims[0] + plot_lims[1]) / 2

        if scale_axis:
            alpha_vec = self.recoverParams(initK, ["prob"])["prob"]
            alpha_vec = alpha_vec[self.order_x]
        else:
            alpha_vec = [1/initK for _ in range(initK)]

        cumul_posn = np.concatenate(([0], np.cumsum(alpha_vec)))
        centers = (cumul_posn[1:] + cumul_posn[:-1]) / 2

        ## Construct matrix
        mat_long = mat_long.assign(
            Z_mod=lambda df: np.where(df['X'] == df['Y'], "--", df['Z']),
            X=lambda df: pd.Categorical(df['X'], categories=self.order_x, ordered=True),
            Y=lambda df: pd.Categorical(df['Y'], categories=self.order_x, ordered=True)
        )

        mat_long['X_idx'] = mat_long['X'].cat.codes
        mat_long['Y_idx'] = mat_long['Y'].cat.codes
        mat_long['X_left'] = cumul_posn[mat_long['X_idx']]
        mat_long['X_right'] = cumul_posn[mat_long['X_idx'] + 1]
        mat_long['X_width'] = mat_long['X_right'] - mat_long['X_left']
        mat_long['X_center'] = (mat_long['X_right'] + mat_long['X_left']) / 2
        mat_long['Y_left'] = cumul_posn[mat_long['Y_idx']]
        mat_long['Y_right'] = cumul_posn[mat_long['Y_idx'] + 1]
        mat_long['Y_width'] = mat_long['Y_right'] - mat_long['Y_left']
        mat_long['Y_center'] = (mat_long['Y_right'] + mat_long['Y_left']) / 2

        ## Visualize
        gridColor = "#000000"
        plt = (p9.ggplot(mat_long, 
                        p9.aes(x='X_center', 
                               y='Y_center', 
                               width='X_width', 
                               height='Y_width', 
                               fill='Z')) + 
            # p9.ggplot(mat_long, p9.aes(x='X', y='Y', fill='Z')) + \
            p9.geom_tile(color=gridColor) + 
            p9.scale_fill_gradient2(
                low="#0000FF",
                midpoint=mid_point,
                high="#FF0000",
                breaks=[plot_lims[0], plot_lims[1]], 
                limits=tuple(plot_lims),
                labels=lambda _: ["Low", "High"],
                name="Similarity"
            )  + \
            p9.theme_bw()  + 
            p9.coord_flip()  + 
            p9.xlab("")  + 
            p9.ylab("")  + 
            # p9.scale_x_discrete(labels=self.order_x)  + \
            # p9.scale_y_discrete(labels=self.order_x)  + \
            p9.scale_x_continuous(
                limits=(0, 1),
                breaks=centers,
                labels=self.order_x,
                expand=(0, 0)
            ) + 
            p9.scale_y_continuous(
                limits=(0, 1),
                breaks=centers,
                labels=self.order_x,
                expand=(0, 0)
            ) + 
            p9.theme(
                # legend_position=legendPosition,
                # axis_text_x=displayAxisFmt,
                # axis_text_y=displayAxisFmt,
                panel_grid=p9.element_blank()
            ))

        return(plt)

    def recoverGroups(self, k):
        if k > self.K:
            raise ValueError("Invalid k provided; {k} > {K}".format(k=k, K=self.K))
        if k == 1:
            return [[x for x in range(self.K)]]
        if k < 1:
            raise ValueError("Invalid k provided; must be > 1")

        ct = [[x] for x in range(self.K)]
        if self.K == k:
            return ct

        for idx in range(self.K-2, k-2, -1):
            res = self.mergeComponents[idx]
            ct[res[0]].extend(ct[res[1]])
            del ct[res[1]]

        return ct

    def recoverParams(self, k, params_to_keep=["prob", "dens", "samplers"]):
        grps = self.recoverGroups(k)
        
        probs = [sum([self.c_prob[y] for y in x]) for x in grps]
        probs = np.array(probs)
        return {
            "prob": probs
        }

    def computePosterior(self, X, k):
        ## Construct group densities
        grps = self.recoverGroups(k)
        merged_clusters = []
        for grp in grps:
            ref_clust = None
            for idx in grp:
                tmp = self.distbn[idx]
                if ref_clust is None:
                    ref_clust = deepcopy(tmp)
                else:
                    ref_clust.merge(tmp)
    
            merged_clusters.append(ref_clust)

        ## Compute Posterior
        return computePosterior(X, merged_clusters)

    def recoverDeltaPmc(self, k):
        grps = self.recoverGroups(k)
        output = np.zeros((k, k))

        ## Consolidate rows
        new_rows = []
        for grp in grps:
            nr = np.sum(self.mat[grp, :], 0)
            nr[grp] = 0
            new_rows.append(nr)

        ## Dimension should be (k x K)
        new_mat = np.array(new_rows)

        ## Consolidate columns
        for idx in range(len(grps)):
            grp = grps[idx]
            nc = np.sum(new_mat[:, grp], 1)
            output[idx, :] = nc
            output[:, idx] = nc

        return output
