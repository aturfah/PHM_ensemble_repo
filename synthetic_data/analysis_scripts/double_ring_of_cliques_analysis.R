# double_ring_of_cliques_analysis.R
#
# The PHM-ensemble pipeline for the double-ring-of-cliques example (36
# cliques / 2 rings). Uses fit_phm_ensemble() from methods_lib.R directly
# -- the same function the six set-of-six examples use -- with its
# default G = 1:25 search range. mclust_auto_extend() (also in
# methods_lib.R, called internally by fit_phm_ensemble()) automatically
# widens that range in steps of G_step=5, up to G_max=100, whenever BIC's
# optimum sits at the boundary of what's been tried -- which is what lets
# every replicate here still reach all 36 cliques starting from G=1:25,
# without this script having to hardcode a wider range itself.
#
# Each replicate's EM is seeded via mclust's own default hc()
# initialization, under the session-wide mclust.options(hcModelName =
# "EII") set near the top of methods_lib.R: mclust's own VVV-seeded
# default gets stuck around K~29-30 for this geometry, while the
# EII-seeded start resolves all 36 cliques cleanly.
#
# M (replicates) and rho (subsample fraction, drawn WITHOUT replacement)
# can be overridden from the command line, mirroring the six set-of-six
# examples:
#   Rscript double_ring_of_cliques_analysis.R --M=5 --rho=0.5
#
# Requires double_ring_of_cliques.rds in ../sim_data/ -- run `make sim`
# (or double_ring_of_cliques_sim.R directly) first if it isn't there yet.
#
# Writes double_ring_of_cliques_phm_result.pdf one directory up (i.e. to
# ~/Dropbox/workspace/phm_ensemble/, alongside the Makefile -- not into
# analysis_scripts/ where this script itself lives), and prints two ARI
# tables (PHM vs. Leiden) comparing both levels of ground-truth
# structure: the fine 36-clique partition (`y`) and the coarse 2-ring
# partition (`ring`).

source("methods_lib.R")
ens_args <- parse_ensemble_args()

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

dat <- readRDS("../sim_data/double_ring_of_cliques.rds")
X <- dat$X
y <- dat$y
ring <- dat$ring
n_blobs <- dat$n_blobs   # 36
n <- nrow(X)
trueK <- n_blobs
label <- "double_ring_of_cliques"

## ================= Stage 1+2: fit_phm_ensemble() -- default G=1:25, auto-extending =================
M <- ens_args$M
rho <- ens_args$rho

cat(sprintf("n=%d, M=%d replicates, rho=%.2f (G starts at 1:25, auto-extends if needed)\n", n, M, rho))
cat(sprintf("START %s\n", format(Sys.time())))

t_start <- Sys.time()
phm <- fit_phm_ensemble(X, M = M, rho = rho, label = label)
total_dt <- as.numeric(Sys.time() - t_start, units = "secs")
cat(sprintf("STAGE1_2_DONE  total_components=%d  total_time=%.1fs\n", length(phm$paramsList), total_dt))

## ================= Stage 3: dendrogram + leaf-ordered merge-height heatmap =================
surface_col <- "white"; ink_col <- "#0b0b0b"
K0 <- length(phm$paramsList)
d  <- length(phm$paramsList[[1]]$mean)
cat(sprintf("K0 (pooled components) = %d, n = %d, d = %d, true K = %d\n", K0, n, d, trueK))

vis <- PHM:::constructVisData(phm, K = K0, scaleHeights = "d*")
xlab_order <- as.integer(vis$xlab)
pos_of <- setNames(seq_along(xlab_order), xlab_order)

sorted_heights <- sort(vis$labels$y)
m <- K0 - trueK
cutYK <- mean(c(sorted_heights[m], sorted_heights[m + 1]))
m2 <- K0 - 2
cutY2 <- mean(c(sorted_heights[m2], sorted_heights[m2 + 1]))
cat(sprintf("cutYK (K=%d) = %.4f, cutY2 (K=2) = %.4f\n", trueK, cutYK, cutY2))

base_theme <- theme_minimal(base_size = 11) +
  theme(plot.background  = element_rect(fill = surface_col, color = NA),
        panel.background = element_rect(fill = surface_col, color = NA),
        panel.grid = element_blank())
xlims <- c(0.5, K0 + 0.5)
common_x <- function() scale_x_continuous(limits = xlims, expand = c(0, 0), breaks = NULL)
common_y <- function() scale_y_reverse(limits = rev(xlims), expand = c(0, 0), breaks = NULL)

p_dendro <- ggplot(vis$df, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment(color = "#3a3a38", linewidth = 0.5) +
  geom_hline(yintercept = cutYK, linetype = "dashed", color = "#c0392b", linewidth = 0.55) +
  geom_hline(yintercept = cutY2, linetype = "dashed", color = "#2a78d6", linewidth = 0.55) +
  common_x() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = NULL, title = "PHM results for double ring of cliques") +
  base_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5, colour = ink_col, margin = margin(b = 8)))

mergeVals <- phm$mergeValues
mergeComps <- phm$mergeComps
label_map <- as.list(seq_len(K0))
merge_matrix <- matrix(NA_real_, K0, K0)
for (xk in K0:2) {
  mc <- mergeComps[[xk]]
  grid <- expand.grid(label_map[[mc[1]]], label_map[[mc[2]]])
  merge_matrix[grid[, 1], grid[, 2]] <- mergeVals[xk]
  merge_matrix[grid[, 2], grid[, 1]] <- mergeVals[xk]
  label_map[[mc[1]]] <- c(label_map[[mc[1]]], label_map[[mc[2]]])
  label_map[[mc[2]]] <- NULL
}
fillScaleFunc <- function(v) -1 * log1p(1 / v)

heat_long <- expand.grid(i = seq_len(K0), j = seq_len(K0))
heat_long$val <- merge_matrix[cbind(heat_long$i, heat_long$j)]
heat_long$Z <- fillScaleFunc(heat_long$val)
heat_long$Z[heat_long$i == heat_long$j] <- NA
heat_long$posx <- pos_of[as.character(heat_long$i)]
heat_long$posy <- pos_of[as.character(heat_long$j)]

plot_lims <- fillScaleFunc(range(heat_long$val, na.rm = TRUE))
p_heat <- ggplot(heat_long, aes(x = posx, y = posy, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        limits = plot_lims, midpoint = mean(plot_lims),
                        na.value = "grey50", guide = "none") +
  common_x() + common_y() +
  base_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

combined <- p_dendro / p_heat + plot_layout(heights = c(3, 6)) &
  theme(plot.background = element_rect(fill = surface_col, color = NA))

out_pdf <- "../double_ring_of_cliques_phm_result.pdf"
ggsave(out_pdf, combined, width = 12, height = 13, bg = surface_col, device = grDevices::cairo_pdf)
cat(sprintf("WROTE %s\n", out_pdf))

## ================= Stage 4: PHM vs. Leiden, at both levels of structure =================
## Fine structure: the 36 individual cliques (`y`). Coarse structure: the
## 2 rings each clique belongs to (`ring`). PHM's point-level partition
## comes from cutting its own merge tree at the target K (phm_partition());
## Leiden's comes from a (k, resolution) grid search targeting that same K
## (run_leiden_grid(), both from methods_lib.R).
cat(sprintf("\n[%s] computing PHM partition + Leiden grid search at K=%d (fine) and K=2 (coarse) ...\n",
            label, trueK))

r_phm_fine     <- phm_partition(phm, X, y, trueK)
r_phm_coarse   <- phm_partition(phm, X, ring, 2)
r_leiden_fine   <- run_leiden_grid(X, y, trueK, label = label)
r_leiden_coarse <- run_leiden_grid(X, ring, 2, label = label)

print_ari_table(sprintf("%s -- fine structure (%d cliques)", label, trueK),
  c("PHM", "Leiden"), c(r_phm_fine$ari, r_leiden_fine$ari))
print_ari_table(sprintf("%s -- coarse structure (2 rings)", label),
  c("PHM", "Leiden"), c(r_phm_coarse$ari, r_leiden_coarse$ari))

cat("ALL_DONE\n")
