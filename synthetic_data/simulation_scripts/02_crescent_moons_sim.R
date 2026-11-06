# 02_crescent_moons_sim.R
#
# Two moons, scaled to n=3000 total (1500 per moon), same generator,
# seed, and noise level as examples/01_two_moons.R.
#
# set.seed(1) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

set.seed(1)

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

# ---- generate data (mirrors sklearn.datasets.make_moons) ----
n <- 1500          # points per moon -> n_total = 3000
noise_sd <- 0.07

t1 <- seq(0, pi, length.out = n)
outer <- cbind(cos(t1), sin(t1))

t2 <- seq(0, pi, length.out = n)
inner <- cbind(1 - cos(t2), 1 - sin(t2) - 0.5)

X <- rbind(outer, inner)
X <- X + matrix(rnorm(length(X), sd = noise_sd), ncol = 2)
y <- c(rep(0, n), rep(1, n))

saveRDS(list(X = X, y = y, n_per_moon = n, noise_sd = noise_sd),
        "two_moons_3000.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
write.table(out_df, "two_moons_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

cat(sprintf("n=%d (per moon=%d)\n", nrow(X), n))

# ---- plot the ground truth ----
png("two_moons_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X, col = point_colors(y), pch = 16, cex = 0.6,
     main = "Two moons (n=3000)", xlab = "x1", ylab = "x2", asp = 1)
dev.off()
