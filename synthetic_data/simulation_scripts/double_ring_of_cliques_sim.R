# double_ring_of_cliques_sim.R
#
# Double ring of cliques (36 blobs: 18 inner + 18 outer), n=3000 total.
# Extracted from ring_of_cliques_2rings_v2.R -- this is just the
# data-generation portion (geometry + simulate + save); the Leiden
# benchmarking/plotting code that followed it in the original script is
# dropped here since it isn't part of generating the simulation data.
#
# set.seed(42) below is the EXACT seed used for the data already
# analyzed -- keep it unchanged to reproduce the same points.
#
# Design goal (v2, fixing v1's flaw): choose the two rings' radii so that
# EVERY blob's own-ring neighbors (both the adjacent one AND the one two
# hops away) are closer than its cross-ring radial partner, for BOTH
# rings -- so "2 rings" is genuinely the coarse-scale ground truth, not
# just a label. Concretely:
#   - inner ring azimuthal edge gap: 3.9 sd
#   - outer ring azimuthal edge gap: whatever falls out at the resulting
#     R_out (same 18 blobs -- comes out ~13 sd here)
#   - radial edge gap (inner blob <-> its outer counterpart): set well
#     ABOVE the outer ring's own gap (edge_gap_sigma_radial = 20), so
#     every blob's nearest and second-nearest neighbors are always on
#     its own ring.
#
# This gives a genuine multi-resolution construction: fine scale = 36
# blobs (cliques), coarse scale = 2 rings.

set.seed(42)

## ================= geometry =================
n_per_ring <- 18
n_blobs    <- 2 * n_per_ring          # 36
sigma      <- 1
edge_radius_sigma <- 3.0
ring_gap_sigma   <- 3.9    # within-ring (inner) azimuthal edge gap
radial_gap_sigma <- 20.0   # cross-ring (radial) edge gap -- deliberately >> any within-ring gap

center_dist_ring   <- 2 * edge_radius_sigma * sigma + ring_gap_sigma * sigma     # 9.9
center_dist_radial <- 2 * edge_radius_sigma * sigma + radial_gap_sigma * sigma   # 26.0

R_in  <- center_dist_ring / (2 * sin(pi / n_per_ring))
R_out <- R_in + center_dist_radial

theta <- 2 * pi * (0:(n_per_ring - 1)) / n_per_ring
centers_in  <- cbind(R_in  * cos(theta), R_in  * sin(theta))
centers_out <- cbind(R_out * cos(theta), R_out * sin(theta))
centers <- rbind(centers_in, centers_out)
ring_of_blob <- c(rep("inner", n_per_ring), rep("outer", n_per_ring))

## ---- verify the ordering we need: same-ring neighbors (1-hop AND 2-hop)
## must be closer than the cross-ring radial partner, for BOTH rings ----
chord <- function(R, hops) 2 * R * sin(hops * pi / n_per_ring)
gap_from_chord <- function(ch) (ch - 2 * edge_radius_sigma * sigma) / sigma

cat("==== geometry check (edge-to-edge gaps, in units of sigma) ====\n")
cat(sprintf("inner ring : 1-hop gap = %.2f sd | 2-hop gap = %.2f sd\n",
            gap_from_chord(chord(R_in, 1)), gap_from_chord(chord(R_in, 2))))
cat(sprintf("outer ring : 1-hop gap = %.2f sd | 2-hop gap = %.2f sd\n",
            gap_from_chord(chord(R_out, 1)), gap_from_chord(chord(R_out, 2))))
cat(sprintf("radial (inner<->outer, same angle) gap = %.2f sd\n",
            gap_from_chord(R_out - R_in)))
cat(sprintf("R_in=%.3f  R_out=%.3f\n", R_in, R_out))

## ================= simulate 3000 points across 36 blobs =================
n_total <- 3000
base_n <- floor(n_total / n_blobs)
rem <- n_total - base_n * n_blobs
n_per_blob <- rep(base_n, n_blobs)
if (rem > 0) n_per_blob[seq_len(rem)] <- n_per_blob[seq_len(rem)] + 1

X <- matrix(NA_real_, n_total, 2)
y <- integer(n_total)
idx <- 1
for (i in seq_len(n_blobs)) {
  ni <- n_per_blob[i]
  X[idx:(idx + ni - 1), ] <- cbind(rnorm(ni, centers[i, 1], sigma), rnorm(ni, centers[i, 2], sigma))
  y[idx:(idx + ni - 1)] <- i - 1
  idx <- idx + ni
}
ring_of_point <- ring_of_blob[y + 1]

cat(sprintf("\n%d blobs (%d inner + %d outer), %d-%d points each (%d total)\n",
            n_blobs, n_per_ring, n_per_ring, min(n_per_blob), max(n_per_blob), n_total))

saveRDS(list(X = X, y = y, ring = ring_of_point, n_blobs = n_blobs, n_per_ring = n_per_ring),
        "double_ring_of_cliques.rds")

# ---- plain tabular export (space-delimited text) ----
out_df <- data.frame(x1 = X[, 1], x2 = X[, 2], blob = y, ring = ring_of_point)
write.table(out_df, "double_ring_of_cliques_data.txt", sep = " ", row.names = FALSE, quote = FALSE)

cat("saved double_ring_of_cliques.rds and double_ring_of_cliques_data.txt\n")
