# methods_lib.R
#
# Shared clustering-method helpers used by the 0X_*_analysis.R scripts in
# this directory. Source this file from each analysis script rather than
# duplicating ~250 lines of parameter-search machinery seven times over.
#
# Methods covered: PHM (bootstrap ensemble), k-means++, Hierarchical/Ward,
# Hierarchical/single-linkage, mclust (fixed K), HDBSCAN, Leiden, OPTICS.
#
# Convention used throughout (matches the rest of this project): any method
# that takes K directly (k-means++, Ward, single-linkage, mclust, and PHM's
# own tree cut) is given the TRUE K. HDBSCAN, Leiden, and OPTICS instead
# pick their own internal parameter (minPts; (k, resolution); (minPts, xi))
# WITHOUT ever looking at point-level labels -- only the target cluster
# COUNT trueK is used to pick the parameter (the widest contiguous range of
# the parameter that yields exactly trueK clusters, or the closest miss if
# nothing hits exactly). The resulting ARI is always computed only
# afterwards, to report the outcome, never to choose the parameter. HDBSCAN
# and OPTICS are scored on ALL points, with noise counted as its own group
# -- the single canonical number for both.

suppressPackageStartupMessages({
  library(mclust); library(dbscan); library(igraph); library(FNN); library(PHM)
})

## mclust's own default hc-initialization model is "VVV" (verified via
## mclust.options("hcModelName")). Overriding it to "EII" here changes the
## hierarchical-clustering partition that seeds every subsequent Mclust()/
## mclustBIC() EM fit in this session (including inside fit_phm_ensemble()'s
## replicate loop, which never passes its own initialization= argument) --
## it does NOT restrict which covariance model the EM is allowed to
## converge to; that's still controlled separately by modelNames. This is a
## session-wide mclust option, not a per-call parameter, so it applies to
## every example sourcing this file.
##
## Empirically (this project): on double_ring_of_cliques, the VVV-seeded
## default gets every replicate stuck around G~28-29 (never finding the
## true 36 cliques), while the EII-seeded start resolves all 36 cleanly on
## every replicate. On compound (M=5, rho=1, equal weighting, same seeds),
## it gave a modest ARI improvement (0.9210 -> 0.9376) without changing the
## eventual BIC-best model (still VVV on every replicate). It has NOT been
## tested on two_moons or two_spirals, whose true clusters are non-convex,
## so this is not confirmed to help uniformly across all six examples.
mclust.options(hcModelName = "EII")

## ================= command-line M / rho overrides =================

#' Parse `--M=<value>` / `--rho=<value>` command-line arguments, as passed
#' by `make analysis M=... RHO=...` (see the top-level Makefile). Falls
#' back to the given defaults -- which match fit_phm_ensemble()'s own
#' defaults (M=5, rho=1.0) -- for any argument not given, including when a
#' script is run directly with no arguments at all, e.g.
#' `Rscript 01_double_rings_analysis.R`.
parse_ensemble_args <- function(M_default = 5, rho_default = 1.0) {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit) == 0) return(default)
    as.numeric(sub(paste0("^--", name, "="), "", hit[1]))
  }
  list(M = get_arg("M", M_default), rho = get_arg("rho", rho_default))
}

## ================= k-means++ =================

kmeanspp_init <- function(X, k, seed) {
  set.seed(seed)
  n <- nrow(X)
  centers <- matrix(NA_real_, k, ncol(X))
  centers[1, ] <- X[sample.int(n, 1), ]
  d2 <- rowSums(sweep(X, 2, centers[1, ], "-")^2)
  for (i in 2:k) {
    probs <- d2 / sum(d2)
    idx <- sample.int(n, 1, prob = probs)
    centers[i, ] <- X[idx, ]
    d2 <- pmin(d2, rowSums(sweep(X, 2, centers[i, ], "-")^2))
  }
  centers
}

#' Hard-assign each row of X to its nearest row of `centers` (squared
#' Euclidean distance), with NO Lloyd's-algorithm refinement afterward.
#' Used by fit_phm_ensemble(..., stage2_init = "random_from_data"): a
#' single assignment step from randomly chosen data points as centers,
#' unlike "kmeans"/"kmeans++" which iterate kmeans() to convergence.
assign_to_nearest_center <- function(X, centers) {
  d2 <- sapply(seq_len(nrow(centers)), function(j) rowSums(sweep(X, 2, centers[j, ], "-")^2))
  max.col(-d2, ties.method = "first")
}

#' Best of n_seeds k-means++ restarts (by lowest total within-cluster SS).
run_kmeanspp <- function(X, y, trueK, n_seeds = 20) {
  best_km <- NULL
  for (seed in 1:n_seeds) {
    init <- kmeanspp_init(X, trueK, seed)
    km <- kmeans(X, centers = init, algorithm = "Hartigan-Wong", iter.max = 100)
    if (is.null(best_km) || km$tot.withinss < best_km$tot.withinss) best_km <- km
  }
  list(cluster = best_km$cluster, ari = adjustedRandIndex(best_km$cluster, y))
}

## ================= hierarchical (Ward / single-linkage) =================

run_ward <- function(X, y, trueK) {
  hc <- hclust(dist(X), method = "ward.D2")
  cl <- cutree(hc, k = trueK)
  list(cluster = cl, ari = adjustedRandIndex(cl, y))
}

run_single_linkage <- function(X, y, trueK) {
  hc <- hclust(dist(X), method = "single")
  cl <- cutree(hc, k = trueK)
  list(cluster = cl, ari = adjustedRandIndex(cl, y))
}

## ================= mclust, fixed K =================

run_mclust_fixed <- function(X, y, trueK) {
  mcl <- Mclust(X, G = trueK, verbose = FALSE)
  list(cluster = mcl$classification, ari = adjustedRandIndex(mcl$classification, y), modelName = mcl$modelName)
}

## ================= HDBSCAN: widest-plateau minPts selection =================

#' @return list(cluster, ari, ari_clustered_only, minPts, K_found, noise_n, noise_pct, scan_df)
run_hdbscan_widest_plateau <- function(X, y, trueK,
                                        minPts_grid = 3:200,
                                        label = "example") {
  n <- nrow(X)
  Ks <- integer(length(minPts_grid))
  noise_ns <- integer(length(minPts_grid))
  cl_cache <- vector("list", length(minPts_grid))
  for (i in seq_along(minPts_grid)) {
    res <- hdbscan(X, minPts = minPts_grid[i])
    cl <- res$cluster
    Ks[i] <- length(unique(cl[cl != 0]))
    noise_ns[i] <- sum(cl == 0)
    cl_cache[[i]] <- cl
  }
  scan_df <- data.frame(minPts = minPts_grid, K_found = Ks, noise_n = noise_ns,
                         noise_pct = 100 * noise_ns / n)

  hit <- which(Ks == trueK)
  if (length(hit) == 0) {
    scan_df$dist_to_trueK <- abs(scan_df$K_found - trueK)
    fallback <- scan_df[order(scan_df$dist_to_trueK, scan_df$minPts), ][1, ]
    cat(sprintf("[%s] HDBSCAN: no minPts on the grid [%d, %d] hit K=%d exactly; nearest was K=%d at minPts=%d\n",
                label, min(minPts_grid), max(minPts_grid), trueK, fallback$K_found, fallback$minPts))
    minPts_final <- fallback$minPts
  } else {
    is_hit <- Ks == trueK
    run_id <- cumsum(c(1, diff(is_hit) != 0))
    runs <- split(seq_along(Ks), run_id)
    best_run <- NULL; best_width <- -1
    for (r in runs) {
      if (!is_hit[r[1]]) next
      width <- length(r)
      if (width > best_width) { best_width <- width; best_run <- r }
    }
    minPts_final <- minPts_grid[best_run[ceiling(length(best_run) / 2)]]
    cat(sprintf("[%s] HDBSCAN: widest K=%d plateau spans minPts [%d, %d] (width=%d) -> picked minPts=%d\n",
                label, trueK, minPts_grid[min(best_run)], minPts_grid[max(best_run)], best_width, minPts_final))
  }

  idx_final <- which(minPts_grid == minPts_final)
  cl_final <- cl_cache[[idx_final]]
  K_found <- length(unique(cl_final[cl_final != 0]))
  noise_n <- sum(cl_final == 0)

  ari <- adjustedRandIndex(cl_final, y)
  keep <- cl_final != 0
  ari_clustered_only <- if (sum(keep) > 1) adjustedRandIndex(cl_final[keep], y[keep]) else NA

  list(cluster = cl_final, ari = ari, ari_clustered_only = ari_clustered_only,
       minPts = minPts_final, K_found = K_found, noise_n = noise_n,
       noise_pct = 100 * noise_n / n, scan_df = scan_df)
}

## ================= Leiden: (k, resolution) grid search =================

build_knn_graph <- function(X, k) {
  nn <- get.knn(X, k = k)
  n <- nrow(X)
  sigma <- nn$nn.dist[, k]
  ii <- rep(seq_len(n), each = k)
  jj <- as.vector(t(nn$nn.index))
  dd <- as.vector(t(nn$nn.dist))
  sim <- exp(-dd^2 / (2 * sigma[ii] * sigma[jj]))
  el <- data.frame(i = ii, j = jj, w = sim)
  key <- ifelse(el$i < el$j, paste(el$i, el$j), paste(el$j, el$i))
  el$key <- key
  el_agg <- aggregate(w ~ key, data = el, FUN = max)
  ij <- do.call(rbind, strsplit(el_agg$key, " "))
  graph_from_data_frame(data.frame(from = ij[, 1], to = ij[, 2], weight = el_agg$w),
                         directed = FALSE, vertices = data.frame(name = as.character(seq_len(n))))
}

#' @return list(cluster, ari, k, resolution, K_found, sweep_df)
run_leiden_grid <- function(X_graph, y, trueK,
                             k_vals = c(10, 15, 20, 30, 50),
                             res_vals = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2,
                                          0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0),
                             seed = 1, label = "example") {
  sweep_rows <- list()
  graphs <- list()
  cl_cache <- list()
  for (k in k_vals) {
    set.seed(seed)
    g <- build_knn_graph(X_graph, k)
    graphs[[as.character(k)]] <- g
    for (res in res_vals) {
      set.seed(seed)
      cl <- cluster_leiden(g, objective_function = "modularity",
                            weights = E(g)$weight, resolution_parameter = res,
                            n_iterations = 10)$membership
      K_found <- length(unique(cl))
      ari <- adjustedRandIndex(cl, y)
      sweep_rows[[length(sweep_rows) + 1]] <- data.frame(k = k, resolution = res, K_found = K_found, ARI = ari)
      cl_cache[[paste(k, res)]] <- cl
    }
  }
  sweep_df <- do.call(rbind, sweep_rows)

  kK <- sweep_df[sweep_df$K_found == trueK, ]
  if (nrow(kK) > 0) {
    best <- kK[which.max(kK$ARI), ]
    cat(sprintf("[%s] Leiden: best among K=%d settings: k=%d, resolution=%.3g -> ARI=%.4f\n",
                label, trueK, best$k, best$resolution, best$ARI))
  } else {
    best <- sweep_df[which.max(sweep_df$ARI), ]
    cat(sprintf("[%s] Leiden: no setting recovered exactly K=%d; using overall-best ARI: k=%d, resolution=%.3g -> K_found=%d, ARI=%.4f\n",
                label, trueK, best$k, best$resolution, best$K_found, best$ARI))
  }

  cl_final <- cl_cache[[paste(best$k, best$resolution)]]
  list(cluster = cl_final, ari = best$ARI, k = best$k, resolution = best$resolution,
       K_found = best$K_found, sweep_df = sweep_df)
}

## ================= OPTICS: widest-plateau (minPts, xi) selection =================

run_optics_widest_plateau <- function(X, y, trueK,
                                       minPts_grid = c(5, 10, 15, 20, 30, 50, 75, 105, 150, 200),
                                       xi_grid = seq(0.01, 0.99, by = 0.01),
                                       label = "example") {
  n <- nrow(X)
  best <- NULL
  scan_rows <- list()
  for (mp in minPts_grid) {
    opt <- optics(X, minPts = mp)
    Ks <- rep(NA_integer_, length(xi_grid))
    noise_ns <- rep(NA_integer_, length(xi_grid))
    cl_cache <- vector("list", length(xi_grid))
    for (i in seq_along(xi_grid)) {
      res <- tryCatch(extractXi(opt, xi = xi_grid[i]), error = function(e) NULL)
      if (is.null(res) || length(res$cluster) != n) next
      cl <- res$cluster
      Ks[i] <- length(unique(cl[cl != 0]))
      noise_ns[i] <- sum(cl == 0)
      cl_cache[[i]] <- cl
    }
    scan_rows[[as.character(mp)]] <- data.frame(minPts = mp, xi = xi_grid, K_found = Ks, noise_n = noise_ns)

    is_hit <- !is.na(Ks) & Ks == trueK
    if (!any(is_hit)) next
    run_id <- cumsum(c(1, diff(is_hit) != 0))
    runs <- split(seq_along(Ks), run_id)
    for (r in runs) {
      if (!is_hit[r[1]]) next
      width <- length(r)
      mid <- r[ceiling(length(r) / 2)]
      noise_pct_mid <- 100 * noise_ns[mid] / n
      better <- is.null(best) ||
        (width > best$width) ||
        (width == best$width && noise_pct_mid < best$noise_pct)
      if (better) {
        best <- list(minPts = mp, width = width, xi = xi_grid[mid],
                     cluster = cl_cache[[mid]], K_found = Ks[mid], noise_pct = noise_pct_mid)
      }
    }
  }
  scan_df <- do.call(rbind, scan_rows)

  if (is.null(best)) {
    valid <- scan_df[!is.na(scan_df$K_found), ]
    valid$dist_to_trueK <- abs(valid$K_found - trueK)
    fb <- valid[order(valid$dist_to_trueK, valid$minPts, valid$xi), ][1, ]
    opt <- optics(X, minPts = fb$minPts)
    cl <- extractXi(opt, xi = fb$xi)$cluster
    cat(sprintf("[%s] OPTICS: no (minPts,xi) hit K=%d exactly; nearest K=%d at minPts=%d, xi=%.3g\n",
                label, trueK, fb$K_found, fb$minPts, fb$xi))
    minPts_final <- fb$minPts; xi_final <- fb$xi
  } else {
    minPts_final <- best$minPts; xi_final <- best$xi; cl <- best$cluster
    cat(sprintf("[%s] OPTICS: widest K=%d plateau at minPts=%d (width=%d) -> xi=%.3g, noise=%.1f%%\n",
                label, trueK, minPts_final, best$width, xi_final, 100 * sum(cl == 0) / n))
  }
  list(cluster = cl, ari = adjustedRandIndex(cl, y), minPts = minPts_final, xi = xi_final)
}

## ================= mclust: BIC search with automatic G-range extension =================

#' Mclust's own (model, G) BIC search, but if the best cell lands on the
#' UPPER boundary of G_init (a sign the true BIC optimum may lie above the
#' range tried), extend the search upward in steps of G_step and merge the
#' new BIC cells onto the ones already computed, via mclust's own
#' mclustBICupdate(). This repeats until the best cell is no longer at the
#' boundary of the accumulated range, or the hard cap G_max is reached.
#'
#' This is NOT an approximation of a wider direct search -- it IS one, just
#' computed incrementally. Verified empirically: Mclust(data, x = merged_bic)
#' is identical to Mclust(data, G = min(G_init):G_hi_final) in every field
#' that matters downstream (G, modelName, loglik, classification, z,
#' parameters, bic, and constructPmcParamsMclust()'s output), because
#' mclustBICupdate() just concatenates BIC matrices computed on the same
#' data and Mclust(x=) only ever reads off the best cell -- no EM is re-run
#' for any G value that was already tried in an earlier round.
#'
#' @param G_init  initial G range to try, e.g. 1:20.
#' @param G_step  how many additional G values to add per extension round.
#' @param G_max   hard cap on G -- extension stops here even if the BIC
#'                optimum is still at the boundary, so a degenerate/
#'                pathological dataset can't trigger an unbounded search.
#' @param label   only used in progress messages when verbose = TRUE.
#' @return an ordinary Mclust object for the best (model, G) found over the
#'         (possibly extended) range -- same fields as a plain Mclust() call.
mclust_auto_extend <- function(data, G_init = 1:25, modelNames = NULL,
                                G_step = 5, G_max = 100,
                                label = "mclust", verbose = TRUE, ...) {
  G_hi <- max(G_init)
  bic  <- mclust::mclustBIC(data, G = G_init, modelNames = modelNames, verbose = FALSE, ...)

  repeat {
    bestG <- as.integer(sub("^.*,([0-9]+)$", "\\1", names(mclust::pickBIC(bic, 1))))
    if (bestG < G_hi || G_hi >= G_max) break   # interior optimum, or hit the safety cap

    new_hi <- min(G_hi + G_step, G_max)
    if (verbose) {
      cat(sprintf("[%s] mclust_auto_extend: BIC-optimal G=%d at upper boundary of G<=%d -- extending search to %d:%d\n",
                  label, bestG, G_hi, G_hi + 1, new_hi))
    }
    bic_ext <- tryCatch(
      mclust::mclustBIC(data, G = (G_hi + 1):new_hi, modelNames = modelNames, verbose = FALSE, ...),
      error = function(e) NULL)
    if (is.null(bic_ext)) {
      if (verbose) {
        cat(sprintf("[%s] mclust_auto_extend: extension to %d:%d failed -- stopping with G<=%d\n",
                    label, G_hi + 1, new_hi, G_hi))
      }
      break
    }
    bic  <- mclust::mclustBICupdate(bic, bic_ext)
    G_hi <- new_hi
  }

  if (verbose && bestG == G_max) {
    cat(sprintf("[%s] mclust_auto_extend: stopped at hard cap G_max=%d with BIC still best at the boundary -- consider raising G_max\n",
                label, G_max))
  }

  mclust::Mclust(data, x = bic)
}

## ================= PHM: ensemble fit from scratch =================

#' Convert an mclust me() EM-fit result (NOT a full Mclust object) into the
#' PHM-compatible per-component params list (prob/mean/var/class).
build_params_from_em <- function(em_fit) {
  K <- length(em_fit$parameters$pro)
  lapply(seq_len(K), function(j) list(
    prob = em_fit$parameters$pro[j],
    mean = em_fit$parameters$mean[, j, drop = FALSE],
    var  = em_fit$parameters$variance$sigma[, , j, drop = FALSE],
    class = j))
}

#' Fit a PHM ensemble of M replicates.
#'
#' Data per replicate: a subsample of size round(rho*n), drawn WITHOUT
#' replacement (NOT a bootstrap resample). At rho == 1 this is a draw of
#' all n rows without replacement, i.e. a full permutation of X -- every
#' row is included exactly once, so it is equivalent to using X as-is;
#' no separate "use all of X" case is needed.
#'
#' (model, G) search (Stage 1, always run): mclust's own default (hc-init)
#' BIC search over G, with the range automatically extended upward (see
#' mclust_auto_extend(), above) whenever the BIC-best G lands on the upper
#' boundary of the range tried for that replicate -- so a too-narrow G
#' upper bound no longer silently truncates the result.
#'
#' Stage 2 (optional, controlled by stage2_init) is where ensemble
#' diversity comes from when it runs: the chosen (model, G) is RE-FIT by
#' EM (mclust::me()) starting from a fresh, replicate-specific random
#' initialization, converted to a responsibility matrix via unmap().
#' Different random starts can converge to different local optima of the
#' same (model, G) -- exactly the mechanism this project's local-optima
#' analysis demonstrated directly (same data, same (model,G), only the EM
#' starting point differs). This is why M > 1 is still meaningful even
#' when rho = 1 (identical data every replicate): the randomness comes
#' from the EM initialization, not from the data. Retries with a fresh
#' random seed if a given start produces a degenerate fit (e.g. an empty
#' component), up to n_restart_tries times, before falling back to the
#' plain Stage-1 fit for that replicate.
#'
#' @param modelNames mclust covariance-model family/families to restrict the
#'   Stage 1 (model, G) BIC search to, e.g. "EII" (spherical, equal volume)
#'   or "VVV" (fully unconstrained) -- or a vector of several, in which case
#'   mclust still searches over G and picks the best model among the ones
#'   listed. NULL (default) imposes no restriction: mclust searches its
#'   full default set of covariance families (14 models when ncol(X) > 1,
#'   or c("E","V") when ncol(X) == 1) and picks whichever wins on BIC. This
#'   is forwarded as-is to mclust_auto_extend()'s mclustBIC()/
#'   mclustBICupdate() calls, so an extended search (see mclust_auto_extend()
#'   above) stays restricted to the same family/families throughout. Whatever
#'   model mclust settles on for a given replicate (mcl0$modelName) is also
#'   what Stage 2's EM re-fit (mclust::me()) uses for that replicate, so
#'   restricting modelNames here fixes the covariance family for both
#'   stages, not just Stage 1.
#' @param stage2_init Stage 2 random-restart method, or NULL to skip Stage
#'   2 entirely (the former fit_phm_ensemble_fast() behavior: each
#'   replicate's contribution is Stage 1's own converged fit, used as-is).
#'   One of:
#'   - NULL           : no Stage 2 (default).
#'   - "kmeans"       : base R kmeans(centers = Gk) -- Gk random data rows
#'                      as starting centers (Forgy init), run to
#'                      convergence. Relies on the ambient RNG stream
#'                      (not explicitly reseeded per attempt), matching
#'                      this project's original fit_phm_ensemble().
#'   - "kmeans++"     : kmeanspp_init() centers, run to convergence via
#'                      kmeans(). Explicitly reseeded each attempt
#'                      (seeds[idx] * 100 + attempt) for reproducibility.
#'   - "random"       : each point independently assigned a uniformly
#'                      random hard label in 1..Gk -- no spatial
#'                      information at all, no clustering run (cheapest
#'                      possible random start; akin to sklearn
#'                      GaussianMixture's init_params = "random", except
#'                      hard-labeled rather than soft responsibilities).
#'                      Explicitly reseeded each attempt.
#'   - "random_from_data" : Gk random data rows used directly as centers,
#'                      with a single nearest-center assignment pass
#'                      (assign_to_nearest_center()) and NO kmeans
#'                      iteration -- spatially informed but cheaper than
#'                      "kmeans"/"kmeans++" (akin to sklearn's
#'                      init_params = "random_from_data"). Explicitly
#'                      reseeded each attempt.
#'
#' @param replicate_weighting How much total probability mass each of the M
#'   replicates contributes to the pooled ensemble (before PHMv2 is fit).
#'   One of:
#'   - "equal"      : default. Every replicate contributes exactly 1/M of
#'                    the total mass, regardless of how many mixture
#'                    components it has or how well it fits -- the
#'                    original behavior (each replicate's own component
#'                    probs, which already sum to 1, are simply rescaled
#'                    by 1/M).
#'   - "likelihood" : each replicate is weighted by its own fit's
#'                    log-likelihood (mcl0$loglik normally, or the
#'                    Stage-2 EM's emr$loglik when Stage 2 ran and
#'                    succeeded for that replicate), via a softmax over
#'                    the M replicate-level log-likelihoods. A replicate
#'                    whose fit explains its (sub)sample better gets
#'                    proportionally more say in the pooled ensemble.
#'                    These are joint log-densities over an entire
#'                    (sub)sample of n rows, so raw log-likelihoods run
#'                    into the thousands-negative range; exponentiating
#'                    one directly underflows to exactly 0 in double
#'                    precision (exp() of anything below roughly -745
#'                    does) long before the *relative* weights would
#'                    differ meaningfully. The fix is the standard
#'                    log-sum-exp shift: subtract the largest
#'                    log-likelihood among the M replicates before
#'                    exponentiating, then normalize so the weights sum
#'                    to 1. This is an exact algebraic identity --
#'                    exp(a)/sum(exp(a)) == exp(a - m)/sum(exp(a - m))
#'                    for any constant m, so shifting by m = max(a)
#'                    changes nothing mathematically, it just keeps
#'                    every exponentiated term in a safe numeric range
#'                    (the largest term becomes exp(0) = 1 instead of a
#'                    number that would otherwise underflow).
#'   - "bic"        : like "likelihood", but weights by the BIC-penalized
#'                    log-likelihood instead of the raw one, so a
#'                    replicate doesn't get more say just because
#'                    mclust_auto_extend() happened to let it settle on
#'                    a larger G (more components almost always means a
#'                    higher raw loglik, whether or not that complexity
#'                    is actually warranted -- BIC weighting corrects
#'                    for that). mclust's own bic() is defined as
#'                    2*loglik - npar*log(n) (so, like loglik, bigger is
#'                    better -- confirmed directly: for a 3-component
#'                    VEI fit mclust::bic(modelName="VEI", loglik=
#'                    m$loglik, n=150, d=2, G=3, equalPro=FALSE) returns
#'                    exactly m$bic, -868.997689, matching
#'                    2*(-404.435033) - 12*log(150)), so the "penalized
#'                    log-likelihood" on the same scale as plain loglik
#'                    is bic/2. The weights are then a softmax over
#'                    these M values, log-sum-exp-shifted exactly as in
#'                    "likelihood" above. Each replicate's BIC is read
#'                    straight off mcl0$bic when Stage 2 didn't run (or
#'                    fell back), since that is already the BIC mclust
#'                    itself computed for the chosen (model, G) cell;
#'                    when Stage 2's EM re-fit succeeded, the BIC is
#'                    recomputed for that fit via mclust::bic(modelName,
#'                    emr$loglik, n, d, Gk, equalPro = FALSE) -- the same
#'                    formula, with equalPro = FALSE throughout since
#'                    this function never changes mclust's default
#'                    there.
#'
#' @return the fitted PHMv2 object (phm$paramsList has K0 pooled components).
fit_phm_ensemble <- function(X, M = 5, rho = 1.0, G = 1:25, modelNames = NULL,
                              replicate_weighting = c("equal", "likelihood", "bic"),
                              seed_base = 1000,
                              numCores = 2, mcSamples = 2e5, batchSize = 1e5,
                              stage2_init = NULL, n_restart_tries = 5,
                              G_step = 5, G_max = 100, extend_verbose = TRUE,
                              label = "example") {
  replicate_weighting <- match.arg(replicate_weighting)
  stage2_methods <- c("kmeans", "kmeans++", "random", "random_from_data")
  if (!is.null(stage2_init) && !stage2_init %in% stage2_methods) {
    stop(sprintf("stage2_init must be NULL or one of: %s (got '%s')",
                  paste(stage2_methods, collapse = ", "), stage2_init))
  }
  if (!is.null(modelNames)) {
    valid_models <- c("E", "V", mclust.options("emModelNames"))
    bad <- setdiff(modelNames, valid_models)
    if (length(bad) > 0) {
      stop(sprintf("modelNames contains unrecognized mclust model name(s): %s (valid: %s)",
                    paste(bad, collapse = ", "), paste(valid_models, collapse = ", ")))
    }
  }

  n <- nrow(X)
  set.seed(seed_base)
  seeds <- sample.int(1e7, M)

  all_params <- vector("list", M)
  replicate_loglik <- numeric(M)
  replicate_bic <- numeric(M)
  for (idx in seq_len(M)) {
    set.seed(seeds[idx])

    ## ---- data for this replicate ----
    ## A single without-replacement draw of round(rho*n) rows. At rho == 1 this
    ## draws all n rows without replacement -- a full permutation of X, i.e.
    ## every row exactly once -- so it needs no special case: row order never
    ## affects any downstream step (Mclust/kmeans/me all treat rows as an
    ## unordered collection), so "all of X, permuted" behaves identically to
    ## "all of X, as-is".
    sub_idx <- sample.int(n, round(rho * n), replace = FALSE)
    sub_dat <- X[sub_idx, ]

    ## ---- Stage 1 (always run): (model, G) chosen via mclust's own
    ## default (hc-init) BIC search, auto-extending the G range upward if
    ## the optimum sits at its boundary (see mclust_auto_extend() above) ----
    mcl0 <- mclust_auto_extend(sub_dat, G_init = G, modelNames = modelNames,
                                G_step = G_step, G_max = G_max,
                                label = paste0(label, " rep", idx), verbose = extend_verbose)
    Gk <- mcl0$G
    modelName <- mcl0$modelName

    if (is.null(stage2_init)) {
      ## ---- no Stage 2: Stage 1's own converged fit, used as-is ----
      params <- constructPmcParamsMclust(mcl0)
      replicate_loglik[idx] <- mcl0$loglik
      replicate_bic[idx] <- mcl0$bic
      cat(sprintf("[%s] PHM replicate %d/%d: G=%d (model=%s), loglik=%.2f (no Stage 2)\n",
                  label, idx, M, Gk, modelName, mcl0$loglik))
    } else {
      ## ---- Stage 2: re-fit EM from a FRESH random initialization at the
      ## chosen (model, G) -- see stage2_init options above. Retry with a
      ## new random seed if a given start produces a degenerate fit, up to
      ## n_restart_tries times, before falling back to the plain mcl0 fit. ----
      emr <- NULL
      for (attempt in seq_len(n_restart_tries)) {
        z0 <- switch(stage2_init,
          "kmeans" = {
            km <- kmeans(sub_dat, centers = Gk, nstart = 1, iter.max = 100)
            unmap(km$cluster, groups = seq_len(Gk))
          },
          "kmeans++" = {
            set.seed(seeds[idx] * 100 + attempt)
            init <- kmeanspp_init(sub_dat, Gk, seeds[idx] * 100 + attempt)
            km <- kmeans(sub_dat, centers = init, iter.max = 100)
            unmap(km$cluster, groups = seq_len(Gk))
          },
          "random" = {
            set.seed(seeds[idx] * 100 + attempt)
            labels <- sample.int(Gk, nrow(sub_dat), replace = TRUE)
            unmap(labels, groups = seq_len(Gk))
          },
          "random_from_data" = {
            set.seed(seeds[idx] * 100 + attempt)
            ctr_idx <- sample.int(nrow(sub_dat), Gk)
            labels <- assign_to_nearest_center(sub_dat, sub_dat[ctr_idx, , drop = FALSE])
            unmap(labels, groups = seq_len(Gk))
          })
        emr <- tryCatch(mclust::me(modelName = modelName, data = sub_dat, z = z0),
                         error = function(e) NULL)
        if (!is.null(emr) && is.finite(emr$loglik)) break
        emr <- NULL
      }
      if (is.null(emr)) {
        cat(sprintf("[%s] PHM replicate %d/%d: random-restart EM (%s) failed after %d tries, falling back to mclust's own fit\n",
                    label, idx, M, stage2_init, n_restart_tries))
        params <- constructPmcParamsMclust(mcl0)
        replicate_loglik[idx] <- mcl0$loglik
        replicate_bic[idx] <- mcl0$bic
      } else {
        params <- build_params_from_em(emr)
        replicate_loglik[idx] <- emr$loglik
        replicate_bic[idx] <- mclust::bic(modelName = modelName, loglik = emr$loglik,
                                           n = nrow(sub_dat), d = ncol(sub_dat), G = Gk,
                                           equalPro = FALSE)
      }
      cat(sprintf("[%s] PHM replicate %d/%d: G=%d (model=%s)%s\n", label, idx, M, Gk, modelName,
                  if (!is.null(emr)) sprintf(", loglik=%.2f [%s]", emr$loglik, stage2_init) else " (fallback fit)"))
    }
    for (g in seq_along(params)) params[[g]]$class <- paste(params[[g]]$class, idx, sep = "_")
    all_params[[idx]] <- params
  }

  ## ---- replicate weights: "equal" (1/M each, the original behavior),
  ## "likelihood" (softmax over the M replicate-level log-likelihoods),
  ## or "bic" (softmax over the M replicate-level BIC/2 values -- the
  ## BIC-penalized log-likelihood, on the same scale as plain loglik).
  ## Both non-equal schemes use the same log-sum-exp shift for numerical
  ## stability -- see @param replicate_weighting above. Each weight
  ## replaces the flat 1/M scale factor applied to that replicate's own
  ## (already-normalized-to-1) component probabilities, so the M
  ## replicates' total pooled mass still sums to 1 under every scheme. ----
  if (replicate_weighting == "equal") {
    rep_w <- rep(1 / M, M)
  } else {
    score <- switch(replicate_weighting,
      "likelihood" = replicate_loglik,
      "bic"        = replicate_bic / 2)
    shift <- max(score)
    rep_w_unnorm <- exp(score - shift)
    rep_w <- rep_w_unnorm / sum(rep_w_unnorm)
    cat(sprintf("[%s] PHM: replicate weights (%s-based): %s\n",
                label, replicate_weighting, paste(sprintf("%.4g", rep_w), collapse = ", ")))
  }

  ens_params <- do.call(c, all_params)
  rep_idx_of <- rep(seq_len(M), times = vapply(all_params, length, integer(1)))
  for (i in seq_along(ens_params)) ens_params[[i]]$prob <- ens_params[[i]]$prob * rep_w[rep_idx_of[i]]

  cat(sprintf("[%s] PHM: fitting PHMv2 on %d pooled components...\n", label, length(ens_params)))
  PHMv2(paramsList = ens_params, scaling = "alpha", mcSamples = mcSamples,
        batchSize = batchSize, numCores = numCores)
}

## ================= PHM: point-level partition from a fit =================

#' Point-level hard cluster assignment (MAP component, then PHM's K-cut) and
#' ARI vs ground truth, for a fitted PHM ensemble -- the same
#' Cholesky-log-density + tree-cut convention used throughout this project.
phm_partition <- function(phm, X, y, trueK) {
  K0 <- length(phm$paramsList)
  n  <- nrow(X); d <- ncol(X)
  log_dens <- matrix(-Inf, n, K0)
  for (k in seq_len(K0)) {
    mu_k <- as.vector(phm$paramsList[[k]]$mean)
    Sigma_k <- matrix(phm$paramsList[[k]]$var, d, d)
    prob_k <- phm$paramsList[[k]]$prob
    ch <- tryCatch(chol(Sigma_k), error = function(e) NULL)
    if (is.null(ch)) next
    logdet <- 2 * sum(log(diag(ch)))
    Xc <- sweep(X, 2, mu_k, "-")
    zz <- t(backsolve(ch, t(Xc), transpose = TRUE))
    maha <- rowSums(zz^2)
    log_dens[, k] <- log(prob_k) - 0.5 * d * log(2 * pi) - 0.5 * logdet - 0.5 * maha
  }
  comp_assign <- max.col(log_dens, ties.method = "first")
  grp <- recoverGroupsv2(phm, trueK)
  point_cluster <- integer(n)
  for (g in seq_along(grp)) point_cluster[comp_assign %in% grp[[g]]] <- g
  list(cluster = point_cluster, ari = adjustedRandIndex(point_cluster, y))
}

## ================= ARI table printing =================

print_ari_table <- function(example_name, method_names, aris) {
  tbl <- data.frame(method = method_names, ARI = round(as.numeric(aris), 4))
  tbl <- tbl[order(-tbl$ARI), ]
  cat(sprintf("\n==== ARI table: %s ====\n", example_name))
  print(tbl, row.names = FALSE)
  invisible(tbl)
}
