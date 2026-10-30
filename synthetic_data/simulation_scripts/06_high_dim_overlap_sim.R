# 06_high_dim_overlap_sim.R
#
# High-dim overlap example, scaled to n=3000 total (750 per cluster,
# 4 clusters), d=60, same seed as examples/12_high_dim_overlap.R, but
# with the cluster centers pulled much closer together (sd 1.4 -> 0.3)
# so the clusters are far less separable, even after PCA denoising.
#
# set.seed(12) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

set.seed(12)

# ---- categorical palette (validated, see dataviz skill reference) ----
palette_cat <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#4a3aa7", "#e34948")
noise_color <- "#898781"
surface_col <- "white"
grid_col    <- "#e1e0d9"
ink_col     <- "#0b0b0b"
point_colors <- function(labels) {
  ifelse(labels < 0, noise_color, palette_cat[(labels %% length(palette_cat)) + 1])
}

# ---- generate data ----
d <- 60
n_per <- 750       # -> 4 * 750 = 3000 total
n_clusters <- 4

center_sd <- 0.5   # was 1.4 -- smaller so clusters overlap more
centers <- matrix(rnorm(n_clusters * d, 0, center_sd), nrow = n_clusters, ncol = d)

X <- matrix(0, nrow = n_clusters * n_per, ncol = d)
y <- integer(n_clusters * n_per)
for (i in seq_len(n_clusters)) {
  rows <- ((i - 1) * n_per + 1):(i * n_per)
  X[rows, ] <- matrix(rnorm(n_per * d, mean = 0, sd = 1), nrow = n_per, ncol = d) +
    matrix(centers[i, ], nrow = n_per, ncol = d, byrow = TRUE)
  y[rows] <- i - 1
}

cat(sprintf("n=%d  d=%d  n_clusters=%d  n_per=%d\n", nrow(X), d, n_clusters, n_per))

# quick illustrative demonstration of distance concentration:
D <- as.matrix(dist(X))
diag(D) <- NA
nn_ratio <- apply(D, 1, function(row) max(row, na.rm = TRUE) / min(row, na.rm = TRUE))
cat(sprintf("mean(farthest/nearest neighbor distance) across points: %.2f\n",
            mean(nn_ratio)))

saveRDS(list(X = X, y = y, d = d, n_per = n_per, n_clusters = n_clusters, centers = centers),
        "high_dim_overlap_3000.rds")

# ---- plain tabular export (space-delimited text): columns x1..x60, y ----
out_df <- as.data.frame(X)
names(out_df) <- paste0("x", seq_len(d))
out_df$y <- y
write.table(out_df, "high_dim_overlap_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

# PCA projection to 2 dimensions, for plotting ONLY
Xc <- scale(X, center = TRUE, scale = FALSE)
pca <- svd(Xc)
X2 <- Xc %*% pca$v[, 1:2]

# ---- plot the ground truth (PCA view) ----
png("high_dim_overlap_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X2, col = point_colors(y), pch = 16, cex = 0.5,
     main = sprintf("High-dim overlap (60-D, n=%d, PCA view)", nrow(X)), xlab = "PC1", ylab = "PC2", asp = 1)
dev.off()

cat("DONE\n")
