# 05_compound_sim.R
#
# "Compound" example (touching clusters of different density + a ring
# enclosing a small blob), scaled to n=3000 total, preserving the
# original 400:400:250:120 ratio as closely as integer counts allow.
# Same seed and distributional parameters as examples/09_compound.R,
# except the blob enclosed within the left-hand ring is now more
# tightly concentrated (sd 0.32 -> 0.15).
#
# set.seed(9) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

set.seed(9)

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

# ---- counts, scaled from 400:400:250:120 (orig total 1170) to sum to 3000 ----
n_dense  <- 1025
n_sparse <- 1026
n_ring   <- 641
n_blob   <- 308
stopifnot(n_dense + n_sparse + n_ring + n_blob == 3000)

# ---- generate data ----
dense <- cbind(rnorm(n_dense, 0, 0.35), rnorm(n_dense, 0, 0.35))
sparse <- cbind(rnorm(n_sparse, 2.3, 0.9), rnorm(n_sparse, 0, 0.9))

theta <- runif(n_ring, 0, 2 * pi)
ring <- cbind(1.15 * cos(theta), 1.15 * sin(theta)) +
  matrix(rnorm(2 * n_ring, 0, 0.06), ncol = 2)
ring <- sweep(ring, 2, c(-4.5, 0), `+`)

blob <- cbind(rnorm(n_blob, -4.5, 0.15), rnorm(n_blob, 0, 0.15))

X <- rbind(dense, sparse, ring, blob)
y <- c(rep(0, n_dense), rep(1, n_sparse), rep(2, n_ring), rep(3, n_blob))

saveRDS(list(X = X, y = y, n_dense = n_dense, n_sparse = n_sparse,
             n_ring = n_ring, n_blob = n_blob),
        "compound_3000.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
write.table(out_df, "compound_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

cat(sprintf("n=%d  (dense=%d, sparse=%d, ring=%d, blob=%d)\n",
            nrow(X), n_dense, n_sparse, n_ring, n_blob))

# ---- plot the ground truth ----
png("compound_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X, col = point_colors(y), pch = 16, cex = 0.6,
     main = "Compound (n=3000)", xlab = "x1", ylab = "x2", asp = 1)
dev.off()
