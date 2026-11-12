require(rhdf5)
require(mclust)
require(PHM)

## Set working directory to be mnist folder location
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grepl("^--file=", args)])
here <- dirname(normalizePath(script_path))
setwd(here)

CORES <- 32

## Ensure files and prepare directories
output_dir <- "output/"
if (!dir.exists(output_dir)) dir.create(output_dir)

density_dir <- file.path(output_dir, "density_replicates")
if (!dir.exists(density_dir)) dir.create(density_dir)

hiha_file <- file.path(output_dir, "human_immune_health_atlas_full.h5ad")
if (!file.exists(hiha_file)) {
    stop("Please place the human_immune_health_atlas_full.h5ad file downloaded from the Human Immune Health Atlas page in the output/ directory")
}

hiha_pcs <- t(h5read(hiha_file, "/obsm/X_pca_harmony"))[, 1:30]

## Run Density Estimation
GMM_REPLICATES <- 50
GMM_N <- 10000
MAX_COMP <- 100
MODELS <- c("VII", "VEI", "EII", "VEI")
gmm_res <- lapply(1:GMM_REPLICATES, function(idx) {
  cat("Replicate", idx, "\n")
  set.seed(idx)
  filename <- paste("gmm_subsamp", "N", GMM_N, "G", MAX_COMP, "repl", idx, "f.RDS", sep="_")
  filename <- file.path(output_dir, filename)

  if (file.exists(filename) & F == T) {
    paramsList <- readRDS(filename)
    return(paramsList)
  } else {
    subsamp_idx <- sample.int(nrow(hiha_pcs), GMM_N, replace=F)
    sub_dat <- hiha_pcs[subsamp_idx, ]
    mcl <- Mclust(sub_dat,
                  G=1:MAX_COMP, 
                  modelNames=MODELS)

    paramsList <- constructPmcParamsMclust(mcl)
    for (g in 1:mcl$G) {
      paramsList[[g]]$class <- paste(paramsList[[g]]$class, idx, sep="_")
    }
    print(mcl$G)
    saveRDS(paramsList, file=filename)
    return(paramsList)
  }
})

params <- do.call(c, gmm_res)
for (idx in 1:length(gmm_res)) {
  gmm_res[[idx]]$prob <- gmm_res[[idx]]$prob / GMM_REPLICATES
}

cat(length(params), "Components\n")

saveRDS(params, file.path(output_dir, "hiha_density.RDS"))