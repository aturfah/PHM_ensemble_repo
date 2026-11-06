## check_packages.R
##
## Verifies every R package used anywhere in this directory -- the
## simulation scripts (01-06 plus double_ring_of_cliques), the shared
## analysis library (methods_lib.R), and double_ring_of_cliques_analysis.R's
## own plotting dependencies -- is installed, and prints install
## instructions for anything that's missing. Invoked by `make check`; not
## meant to be run standalone, but it's plain Rscript so nothing stops you.
##
## Exit status: 0 if everything required is installed, 1 if anything is
## missing (so `make check` fails loudly rather than silently).

cran_pkgs <- c(
  "RANN",    # 04_two_spirals_sim.R only -- a nearest-neighbor distance
             # diagnostic, not used to generate the data itself
  "mclust",  # methods_lib.R: density estimation backing fit_phm_ensemble()
             # and run_mclust_fixed()
  "dbscan",  # methods_lib.R: run_hdbscan_widest_plateau(),
             # run_optics_widest_plateau()
  "igraph",     # methods_lib.R: run_leiden_grid()
  "FNN",        # methods_lib.R: build_knn_graph()
  "ggplot2",    # double_ring_of_cliques_analysis.R: dendrogram + merge-height heatmap
  "patchwork"   # double_ring_of_cliques_analysis.R: combines the two ggplot panels
)

github_pkgs <- list(
  PHM = "aturfah/PHM"   # methods_lib.R: fit_phm_ensemble(), phm_partition() --
                        # the Distinguishability Criterion / PHM algorithms,
                        # not on CRAN
)

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

missing_cran   <- Filter(Negate(is_installed), cran_pkgs)
missing_github <- Filter(function(pkg) !is_installed(pkg), names(github_pkgs))

cat(sprintf("R found: %s\n", R.version.string))
cat(sprintf("Checked %d CRAN package(s) and %d GitHub package(s) ",
            length(cran_pkgs), length(github_pkgs)))
cat(sprintf("required by the simulation and analysis scripts in this directory.\n\n"))

if (length(missing_cran) == 0 && length(missing_github) == 0) {
  cat("All required packages are installed. You're set to run any of the\n")
  cat("0X_*_sim.R / 0X_*_analysis.R scripts, or double_ring_of_cliques_sim.R /\n")
  cat("double_ring_of_cliques_analysis.R, directly.\n")
  quit(status = 0, save = "no")
}

cat("MISSING PACKAGES\n")
cat("=================\n\n")

install_lines <- character(0)

if (length(missing_cran) > 0) {
  cat("CRAN:\n")
  for (p in missing_cran) cat(sprintf("  - %s\n", p))
  cran_call <- sprintf('install.packages(c(%s))',
                        paste(sprintf('"%s"', missing_cran), collapse = ", "))
  cat(sprintf("\nInstall with:\n\n  %s\n\n", cran_call))
  install_lines <- c(install_lines, cran_call)
}

if (length(missing_github) > 0) {
  cat("GitHub (not on CRAN):\n")
  for (p in missing_github) cat(sprintf("  - %s  (from %s)\n", p, github_pkgs[[p]]))
  cat("\nInstall with (requires the 'remotes' package to fetch from GitHub):\n\n")
  if (!is_installed("remotes")) {
    cat('  install.packages("remotes")\n')
    install_lines <- c(install_lines, 'install.packages("remotes")')
  }
  for (p in missing_github) {
    call <- sprintf('remotes::install_github("%s")', github_pkgs[[p]])
    cat(sprintf("  %s\n", call))
    install_lines <- c(install_lines, call)
  }
  cat("\n")
}

cat("All of the above, as a single shell command:\n\n")
cat(sprintf("  Rscript -e '%s'\n\n", paste(install_lines, collapse = "; ")))
cat("Then re-run `make check` to confirm.\n")

quit(status = 1, save = "no")
