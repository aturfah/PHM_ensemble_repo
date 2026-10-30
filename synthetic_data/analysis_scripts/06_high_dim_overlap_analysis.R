# 06_high_dim_overlap_analysis.R
#
# Runs all 8 methods on the high-dim-overlap example and prints the
# resulting ARI table to screen. Requires high_dim_overlap_3000.rds in
# ../sim_data/ -- run 06_high_dim_overlap_sim.R first if it isn't
# there yet.
#
# IMPORTANT: the raw data is 60-dimensional. Per this project's established
# convention (distance concentration bites density-based methods already
# at "modest" dimensionality, well before d=60), EVERY method here --
# including PHM -- clusters on the first 4 PCA scores, not the raw 60-dim
# data. This is the one example where the analysis script does a PCA step
# before calling into methods_lib.R.

source("methods_lib.R")
ens_args <- parse_ensemble_args()

dat <- readRDS("../sim_data/high_dim_overlap_3000.rds")
X_raw <- dat$X; y <- dat$y
trueK <- length(unique(y))
label <- "high_dim_overlap"

## ---- reduce to the first 4 PCs; ALL methods cluster on this, not X_raw ----
Xc <- scale(X_raw, center = TRUE, scale = FALSE)
pca <- svd(Xc)
X <- Xc %*% pca$v[, 1:4]
cat(sprintf("\n===== High-Dim Overlap (n=%d, d_raw=%d -> 4 PCs, true K=%d, M=%d, rho=%.2f) =====\n",
            nrow(X), ncol(X_raw), trueK, ens_args$M, ens_args$rho))

phm <- fit_phm_ensemble(X, M = ens_args$M, rho = ens_args$rho, label = label)
r_phm <- phm_partition(phm, X, y, trueK)

r_km <- run_kmeanspp(X, y, trueK)
r_wd <- run_ward(X, y, trueK)
r_sl <- run_single_linkage(X, y, trueK)
r_mc <- run_mclust_fixed(X, y, trueK)
r_hdb <- run_hdbscan_widest_plateau(X, y, trueK, label = label)
r_lei <- run_leiden_grid(X, y, trueK, label = label)
r_opt <- run_optics_widest_plateau(X, y, trueK, label = label)

print_ari_table("High-Dim Overlap",
  c("PHM", "k-means++", "Hierarchical (Ward)", "Hierarchical (single)", "mclust", "HDBSCAN", "Leiden", "OPTICS"),
  c(r_phm$ari, r_km$ari, r_wd$ari, r_sl$ari, r_mc$ari, r_hdb$ari, r_lei$ari, r_opt$ari))
