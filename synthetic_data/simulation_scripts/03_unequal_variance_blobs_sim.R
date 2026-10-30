# 03_unequal_variance_blobs_sim.R
#
# "Unequal-variance blobs" example (three Gaussian blobs with very
# different spreads: sd = 0.4, 1.2, 2.5), n=3000 total, 1000 per blob
# (exact 3-way split of 3000, same equal-count ratio as the original
# examples/04_unequal_variance_blobs.R at n=900). Same seed, centers,
# and spreads as that original script.
#
# set.seed(4) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

set.seed(4)

palette_cat <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#4a3aa7", "#e34948")
surface_col <- "white"
grid_col    <- "#e1e0d9"
ink_col     <- "#0b0b0b"
point_colors <- function(labels) palette_cat[(labels %% length(palette_cat)) + 1]

n_per <- 1000
centers <- rbind(c(-4, 2), c(3, 3), c(0, -4))
stds <- c(0.4, 1.2, 2.5)

blob <- function(center, n, sd) {
  cbind(rnorm(n, center[1], sd), rnorm(n, center[2], sd))
}

X <- rbind(blob(centers[1, ], n_per, stds[1]),
           blob(centers[2, ], n_per, stds[2]),
           blob(centers[3, ], n_per, stds[3]))
y <- rep(0:2, each = n_per)

saveRDS(list(X = X, y = y, n_per = n_per, centers = centers, stds = stds),
        "unequal_variance_blobs_3000.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
write.table(out_df, "unequal_variance_blobs_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

png("unequal_variance_blobs_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X, col = point_colors(y), pch = 16, cex = 0.5,
     main = sprintf("Unequal-variance blobs, n=%d", nrow(X)), xlab = "x1", ylab = "x2", asp = 1)
dev.off()

cat(sprintf("n=%d total (n_per=%d each; sd = %s)\n", nrow(X), n_per, paste(stds, collapse = ", ")))
cat("DONE\n")
