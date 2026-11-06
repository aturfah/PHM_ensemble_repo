# 01_double_rings_analysis.R
#
# Runs all 8 methods on the double-rings (concentric annuli) example and
# prints the resulting ARI table to screen. Requires
# concentric_annuli_3000.rds in ../sim_data/ -- run
# 01_double_rings_sim.R first if it isn't there yet.

source("methods_lib.R")
ens_args <- parse_ensemble_args()

dat <- readRDS("../sim_data/concentric_annuli_3000.rds")
X <- dat$X; y <- dat$y
trueK <- length(unique(y))
label <- "double_rings"
cat(sprintf("\n===== Double Rings (n=%d, true K=%d, M=%d, rho=%.2f) =====\n",
            nrow(X), trueK, ens_args$M, ens_args$rho))

phm <- fit_phm_ensemble(X, M = ens_args$M, rho = ens_args$rho, label = label)
r_phm <- phm_partition(phm, X, y, trueK)

r_km <- run_kmeanspp(X, y, trueK)
r_wd <- run_ward(X, y, trueK)
r_sl <- run_single_linkage(X, y, trueK)
r_mc <- run_mclust_fixed(X, y, trueK)
r_hdb <- run_hdbscan_widest_plateau(X, y, trueK, label = label)
r_lei <- run_leiden_grid(X, y, trueK, label = label)
r_opt <- run_optics_widest_plateau(X, y, trueK, label = label)

print_ari_table("Double Rings",
  c("PHM", "k-means++", "Hierarchical (Ward)", "Hierarchical (single)", "mclust", "HDBSCAN", "Leiden", "OPTICS"),
  c(r_phm$ari, r_km$ari, r_wd$ari, r_sl$ari, r_mc$ari, r_hdb$ari, r_lei$ari, r_opt$ari))
