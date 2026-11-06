import matplotlib
matplotlib.use('Agg')

import dill
import numpy as np
import anndata as ad
import plotnine as p9
import pandas as pd

from config import output_dir, image_dir, phm_fname, hiha_fname
from helpers import fuzzy_contingency, grid_figure_plot

## Load PHM run
print("==== Loading HIHA PHM ====\n")
if not phm_fname.exists():
    raise RuntimeError("PHM file not found. Run profile_phm before this file!")

with phm_fname.open("rb") as sf:
    phm = dill.load(sf)

VIS_K = 250
# VIS_SUBSAMP = 5e5
VIS_SUBSAMP = 1e4


######################################################
### PHM Dendrogram and Heatmap for overall merging ###
######################################################
print("==== Constructing Cell-type Visualizations ====")

print("\tConstructing overall dendrogram/heatmap")
dend = (phm.dendrogram(initK=VIS_K, scale_heights="d*") +
       p9.theme(figure_size=[8, 4],
                axis_text_x=p9.element_blank(),
                axis_ticks=p9.element_blank())
       )
htmp = (phm.heatmap(initK=VIS_K, scale_heights="d*", scale_axis=True,
                    grid_color="none",
                    grid_size=0) +
       p9.theme(figure_size=[8, 4],
                axis_text_x=p9.element_blank(),
                axis_text_y=p9.element_blank(),
                axis_ticks=p9.element_blank(),
                panel_grid=p9.element_blank()
               ))

dendhtmp_all = (dend | htmp)
dendhtmp_all.save(filename=image_dir.joinpath("hiha_all_dendhtmp.png"),
                 width=8, height=4)

###########################################################
### Posterior matrix for cell-type-specific dendrograms ###
###########################################################
print("\tConstructing posterior matrix for cell type visualization")
adata = ad.read_h5ad(str(hiha_fname), backed="r")
hiha_pca = adata.obsm["X_pca_harmony"][:, :30]

subsamp_file = output_dir.joinpath("vis_ALL_subsamp_idx.npy")
if not subsamp_file.exists():
    np.random.seed(20251104)
    subsamp = np.random.randint(hiha_pca.shape[0], size=int(VIS_SUBSAMP))
    np.save(subsamp_file, subsamp)
else:
    subsamp = np.load(subsamp_file)

posterior_file = output_dir.joinpath("vis_ALL_posterior_{}.npy".format(VIS_K))
if not posterior_file.exists():
    subsamp_post = phm.computePosterior(hiha_pca[subsamp, :], VIS_K)
    np.save(posterior_file, subsamp_post)
else:
    subsamp_post = np.load(posterior_file)


############################################################
### Generate the class-specific dendrograms at cell-type ###
############################################################
print("\tConstructing cell subtype dendrograms")
cell_labels1 = np.asarray(adata.obs["AIFI_L1"])
cell_type_plots = []
subsamp_labels = cell_labels1[subsamp]
fc_table = fuzzy_contingency(subsamp_labels, subsamp_post)
for col_name, col_values in fc_table.items():
    print("\t\tCell Type: {}".format(col_name))
    col_values = list(col_values.to_numpy())
    ctype_dend = (phm.dendrogram(initK=VIS_K,
                                scale_heights="d*",
                                groupColorMax="#FF0000",
                                group_probs=col_values) +
                  p9.ggtitle(col_name) + 
                  p9.theme(axis_text_x=p9.element_blank(),
                           axis_ticks_x=p9.element_blank(),
                           figure_size=[20, 8]))
    cell_type_plots.append(ctype_dend)

subtrees_all = grid_figure_plot(cell_type_plots, 5) & p9.theme(figure_size=[8, 2], title=p9.element_text(size=8))
subtrees_all.save(filename=image_dir.joinpath("hiha_all_celltype_dend.png"),
                 width=8, height=2)

print("\tComplete!\n")

####################################################
### Posterior matrix for T + NK cell dendrograms ###
####################################################
print("==== Constructing T + NK Cell Visualizations ====")

## Retrieve the indices for cells corresponding to these types
target_cells = ["T cell", "NK cell"]
cell_labels2 = np.asarray(adata.obs["AIFI_L2"])

indices = {}
types = {}
for cell in target_cells:
    idx = np.flatnonzero(cell_labels1 == cell)
    indices[cell] = idx
    types[cell] = np.unique(cell_labels2[idx])

subset_idx = np.concatenate([indices[cell] for cell in target_cells])

## Construct the subsample and posterior matrix for these cells
print("\tConstructing posterior matrix for cell type visualization")
subsamp2_file = output_dir.joinpath("vis_TNK_subsamp_idx.npy")
if not subsamp2_file.exists():
    subsamp2 = subset_idx[np.random.randint(len(subset_idx), size=int(VIS_SUBSAMP))]
    np.save(subsamp2_file, subsamp2)
else:
    subsamp2 = np.load(subsamp2_file)

posterior2_file = output_dir.joinpath("vis_TNK_posterior_{}.npy".format(VIS_K))
if not posterior2_file.exists():
    subsamp_post2 = phm.computePosterior(hiha_pca[subsamp2, :], VIS_K)
    np.save(posterior2_file, subsamp_post2)
else:
    subsamp_post2 = np.load(posterior2_file)

subsamp_labels2 = cell_labels2[subsamp2]
fc_table2 = fuzzy_contingency(subsamp_labels2, subsamp_post2)

######################################################################
### Generate class-specific dendrograms for NK and T cell subtypes ###
######################################################################
types["T + NK cell"] = []
types["T + NK cell"].extend(types["T cell"])
types["T + NK cell"].extend(types["NK cell"])

column_indices = {}
for cell in types:
    for t in types[cell]:
        column_indices[t] = list(fc_table2.columns).index(t)

ctype = "T + NK cell"
figures = {
    "ctype": ctype,
    "dend": [],
    "total_dend": [],
    "htmp": []
}
htmp_idx = []
total_weight = None
print("\tConstructing cell subtype dendrograms")
for subctype in types["T + NK cell"]:
    subctype_idx = column_indices[subctype]
    print("\t\tCell Type:", ctype, "|", subctype, "|", subctype_idx)
    col_values = np.asarray(fc_table2.iloc[:, subctype_idx])
    figures["dend"].append(
        (phm.dendrogram(initK=VIS_K,
                       scale_heights="d*",
                       groupColorMax="#FF0000",
                       group_probs=list(col_values)) + 
         p9.ggtitle(subctype) + 
         p9.theme(axis_text_x=p9.element_blank(),
                  axis_ticks=p9.element_blank())
        )
    )
    htmp_idx.extend(np.where(col_values > 0.5)[0])
    if total_weight is None:
        total_weight = col_values
    else:
        total_weight += col_values


print("\tConstructing T+NK dendrogram/heatmap")
figures["htmp"] = phm.heatmap(initK=VIS_K, 
                               scale_heights="d*", 
                               scale_axis=True, 
                               grid_color="none",
                               grid_size=0,
                               subset_index=htmp_idx
                              ) + \
       p9.theme(figure_size=[20, 9],
            axis_text_x=p9.element_blank(),
            axis_text_y=p9.element_blank(),
            axis_ticks=p9.element_blank(),
            panel_grid=p9.element_blank()
           )
figures["total_dend"] = phm.dendrogram(initK=VIS_K,
                       scale_heights="d*",
                       groupColorMax="#FF0000",
                       group_probs=list(total_weight)) + \
        p9.ggtitle(ctype) + \
        p9.theme(axis_text_x=p9.element_blank(),
                 axis_ticks=p9.element_blank())

subtrees_t_nk = grid_figure_plot(figures["dend"], 7) & \
    p9.theme(title=p9.element_text(size=8), figure_size=[8, 2.5])
dendhtmp_t_nk = (figures["total_dend"] | figures["htmp"]) & \
    p9.theme(title=p9.element_text(size=8), figure_size=[8, 3.5])
plt_t_nk = (dendhtmp_t_nk / subtrees_t_nk)  & p9.theme(figure_size=[8, 6])

subtrees_t_nk.save(filename=image_dir.joinpath("hiha_tnk_subtype_dend.png"))
dendhtmp_t_nk.save(filename=image_dir.joinpath("hiha_tnk_dendhtmp.png"))
plt_t_nk.save(filename=image_dir.joinpath("hiha_tnk_combined.png"))
print("\tComplete!\n")
