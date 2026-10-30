# 01_double_rings_sim.R
#
# Two concentric annuli (not thin rings): points uniformly distributed
# over the 2D area of each annulus.
#   inner annulus: radius in [4, 5]
#   outer annulus: radius in [8, 9]
# Total n = 3000, split evenly (1500 / 1500).
#
# set.seed(2) below is the EXACT seed used for the data already analyzed --
# keep it unchanged to reproduce the same points.

set.seed(2)

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
# Uniform sampling over an annulus [r_min, r_max]: radius must be drawn as
# sqrt(runif(n, r_min^2, r_max^2)) so that area (not radius) is uniform.
annulus <- function(r_min, r_max, n) {
  theta  <- runif(n, 0, 2 * pi)
  radius <- sqrt(runif(n, r_min^2, r_max^2))
  cbind(radius * cos(theta), radius * sin(theta))
}

n_total <- 3000
n_inner <- n_total / 2
n_outer <- n_total / 2

r1 <- annulus(4, 5, n_inner)   # inner annulus: radius in [4, 5]
r2 <- annulus(8, 9, n_outer)   # outer annulus: radius in [8, 9]

X <- rbind(r1, r2)
y <- c(rep(0, nrow(r1)), rep(1, nrow(r2)))

saveRDS(list(X = X, y = y, n_inner = n_inner, n_outer = n_outer,
             inner_range = c(4, 5), outer_range = c(8, 9)),
        "concentric_annuli_3000.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
write.table(out_df, "concentric_annuli_3000.txt", sep = " ", row.names = FALSE, quote = FALSE)

cat(sprintf("n=%d  (inner=%d, outer=%d)\n", nrow(X), n_inner, n_outer))
cat(sprintf("inner radius range realized: [%.3f, %.3f]\n", min(sqrt(rowSums(r1^2))), max(sqrt(rowSums(r1^2)))))
cat(sprintf("outer radius range realized: [%.3f, %.3f]\n", min(sqrt(rowSums(r2^2))), max(sqrt(rowSums(r2^2)))))
cat(sprintf("gap between annuli: %.3f\n", min(sqrt(rowSums(r2^2))) - max(sqrt(rowSums(r1^2)))))

# ---- plot the ground truth ----
png("concentric_annuli_3000.png", width = 1600, height = 1400, res = 220, bg = surface_col)
par(bg = surface_col, col.axis = ink_col, col.lab = ink_col, col.main = ink_col,
    fg = grid_col, mar = c(4, 4, 3, 1))
plot(X, col = point_colors(y), pch = 16, cex = 0.6,
     main = "Concentric annuli (n=3000)", xlab = "x1", ylab = "x2", asp = 1)
dev.off()
