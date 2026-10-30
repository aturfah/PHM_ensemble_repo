# 04_two_spirals_sim.R
#
# Two spirals, 1.5 turns each, n=3000 total (1500/arm). Re-tuned so that
# the gap between the two starting points (2*r0) is 1.5x wider than the
# band spacing between the interleaved strands away from the center
# (which, for this r = r0 + k*theta parametrization, equals exactly
# k*pi regardless of r0 or turns -- verified numerically with a huge r0
# so the start point never binds).
#
# k kept at 6.0 (unchanged from the previous version).
#   band  = k*pi        = 18.8496
#   gap   = 1.5 * band   = 28.2743  ->  r0 = gap / 2 = 14.1372
#
# set.seed(7) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

suppressPackageStartupMessages(library(RANN))
set.seed(7)

palette_cat <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#4a3aa7", "#e34948")
surface_col <- "white"
grid_col    <- "#e1e0d9"
ink_col     <- "#0b0b0b"

n <- 1500   # per arm -> 3,000 total
turns <- 1.5
noise_sd <- 0.12
k  <- 6.0
band <- k * pi
gap_at_start <- 1.5 * band
r0 <- gap_at_start / 2

theta <- sqrt(runif(n, 0, 1)) * turns * 2 * pi
r <- r0 + k * theta

x1 <- r * cos(theta) + rnorm(n, 0, noise_sd)
y1 <- r * sin(theta) + rnorm(n, 0, noise_sd)
x2 <- -r * cos(theta) + rnorm(n, 0, noise_sd)
y2 <- -r * sin(theta) + rnorm(n, 0, noise_sd)

X <- rbind(cbind(x1, y1), cbind(x2, y2))
y <- c(rep(0, n), rep(1, n))

saveRDS(list(X = X, y = y, r0 = r0, k = k, n_per_arm = n, turns = turns, noise_sd = noise_sd,
             band_spacing = band, gap_at_start = gap_at_start),
        "two_spirals_3000.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
write.table(out_df, "two_spirals_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

A <- X[y == 0, , drop = FALSE]; B <- X[y == 1, , drop = FALSE]
nn <- nn2(B, A, k = 1)
minDist <- min(nn$nn.dists)
cat(sprintf("r0=%.4f  k=%.2f  turns=%.2f  noise_sd=%.2f\n", r0, k, turns, noise_sd))
cat(sprintf("nominal gap at start (2*r0)     = %.4f\n", 2 * r0))
cat(sprintf("nominal band spacing (k*pi)     = %.4f\n", band))
cat(sprintf("ratio (gap/band)                = %.4f\n", (2 * r0) / band))
cat(sprintf("TRUE realized global min separation (actual sampled points) = %.4f\n", minDist))

point_colors <- function(labels) ifelse(labels == 0, palette_cat[1], palette_cat[2])

png("two_spirals_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X, col = point_colors(y), pch = 16, cex = 0.6,
     main = sprintf("Two spirals (1.5 turns), n=%d, start gap 1.5x band", nrow(X)),
     xlab = "x1", ylab = "x2", asp = 1)
dev.off()

cat("DONE\n")
