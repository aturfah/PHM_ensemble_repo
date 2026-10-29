import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import plotnine as p9
import dill
from pathlib import Path

from sklearn.metrics import adjusted_rand_score
from config import save_dir_global_ensemble, images_directory
from prepare_data import load_tsne, load_raw
from vis_helpers import produce_grid_plot, map_clusters_to_label_mask


proj_tsne = load_tsne()
raw_mnist, sampleY_t = load_raw()

print("\n==== Loading First PHM Ensemble (Global) Replicate ====\n")
reference_file = save_dir_global_ensemble.joinpath("PHM_mclustSubaggALL2_r1.pkl")
if not reference_file.exists():
    raise RuntimeError("File not found: {}. Make sure run_cluster_performance.py has been ran before calling this file!".format(reference_file))

with reference_file.open("rb") as sf:
    phm = dill.load(sf)


#################################
### Dendrogram + Heatmap Plot ###
#################################
print("==== Constructing Dendrogram + Heatmap ====")
TREE_CUTOFF = 20
post_mat = phm.computePosterior(proj_tsne, TREE_CUTOFF)
unique_classes = np.unique(sampleY_t)
class_mask = (sampleY_t[:, None] == unique_classes).astype(float)
M = (class_mask.T @ post_mat) / post_mat.sum(axis=0)
print("\tCompleted Posterior")

## Constructs the Class colored dendrograms
class_plots = []
for class_idx in range(len(unique_classes)):
    class_ = unique_classes[class_idx]
    print("\tDigit: {}".format(class_), end="\r")
    class_prob_vec = list(M[class_idx, :])
    class_plots.append(phm.dendrogram(initK=M.shape[1], 
                                      scale_heights="d*",
                                      group_probs=class_prob_vec, 
                                      groupColorMax="#FF0000") + p9.ggtitle("Digit: {}".format(class_)) +
                      p9.theme(figure_size=[20, 9],
                               legend_position="none",
                               axis_text_x=p9.element_blank()))
print("\tCompleted Class-specific Dendrograms")

## Heatmap showing the overal merging
htmp = phm.heatmap(initK=TREE_CUTOFF, scale_axis=True, grid_size=0, scale_heights="d*",) + p9.theme(
    figure_size=[20, 9], 
    axis_text_x = p9.element_blank(), 
    axis_text_y = p9.element_blank())

## Generally the digits should group (0, 6, 3, 5, 8, 2) VS (4, 9, 7, 1) but the exact ordering does not always correpond
# htmp_txt = htmp + p9.geom_text(
#     mapping=p9.aes(x="x", y="y", label="label", width=None, height=None),
#     data=pd.DataFrame.from_dict([
#         {"x":0.05, "y":0.05, "label": "Digit: 0", "Z": 0},
#         {"x":0.15, "y":0.15, "label": "Digit: 6", "Z": 0},
#         {"x":0.25, "y":0.25, "label": "Digit: 3", "Z": 0},
#         {"x":0.35, "y":0.35, "label": "Digit: 5", "Z": 0},
#         {"x":0.44, "y":0.44, "label": "Digit: 8", "Z": 0},
#         {"x":0.53, "y":0.53, "label": "Digit: 2", "Z": 0},
#         {"x":0.62, "y":0.62, "label": "Digit: 4", "Z": 0},
#         {"x":0.73, "y":0.73, "label": "Digit: 9", "Z": 0},
#         {"x":0.83, "y":0.83, "label": "Digit: 7", "Z": 0},
#         {"x":0.94, "y":0.94, "label": "Digit: 1", "Z": 0}
#     ]),
#     size=16
# )
htmp_txt = htmp
print("\tCompleted Heatmap")

plt_digit_dend = (
    (class_plots[0] | class_plots[6] | class_plots[3] | class_plots[5] | class_plots[8]) / 
    (class_plots[2] | class_plots[4] | class_plots[9] | class_plots[7] | class_plots[1])
) | htmp_txt
print("\tCompleted Composite Plot")

plt_digit_dend.save(
    images_directory.joinpath("mnist_phm_heatmap.png"),
    width=8, height=6)

print("==== Completed Dendrogram + Heatmap ====\n")

#################################
### Dendrogram + Heatmap Plot ###
#################################
print("==== Generating Digit Grid ====")
GRID_WIDTH = 5

## Get digits corresponding to 1, 2, 9, 4, 7. Generally assume each digit has 2 classes as this is what happens
masks_1 = map_clusters_to_label_mask(post_mat, sampleY_t, 1)
masks_2 = map_clusters_to_label_mask(post_mat, sampleY_t, 2)
masks_9 = map_clusters_to_label_mask(post_mat, sampleY_t, 9)
masks_4 = map_clusters_to_label_mask(post_mat, sampleY_t, 4)
masks_7 = map_clusters_to_label_mask(post_mat, sampleY_t, 7)

grid_combined = produce_grid_plot(
    raw_mnist,
    [
        masks_1[0], masks_2[0], masks_9[0], masks_4[0], masks_7[0], 
        masks_1[1], masks_2[1], masks_9[1], masks_4[1], masks_7[1]
    ],
    [
        "Cluster 1a", "Cluster 2a", "Cluster 9a", "Cluster 4a", "Cluster 7a",
        "Cluster 1b", "Cluster 2b", "Cluster 9b", "Cluster 4b", "Cluster 7b"
    ],
    GRID_WIDTH,
    seed=0,
    ncol=5, nrow=2
)  + p9.theme(figure_size=[7, 4.5])


grid_combined.save(
    images_directory.joinpath("mnist_digits.png"),
     width=7, height=4.5)