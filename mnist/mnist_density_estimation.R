require(mclust)
require(PHM)
require(aricode)
require(reticulate)
require(parallel)
require(clue)
require(fastcluster)

rm(list=ls())

np = reticulate::import("numpy")

## Set working directory to be mnist folder location
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grepl("^--file=", args)])
here <- dirname(normalizePath(script_path))
setwd(here)

## Define the locations
ensemble_directory = "output/70k_ensemble_global/"
single_directory = "output/70k_single_PHM"
mclust_directory = "output/70k_mclustG10/"
gt_file = "output/gt_deltapmc_tsne_70k.RData"
fname = "output/MNIST_tsne_70k.npy"
X = np$load(fname)
num_pcs = ncol(X)
labels = np$load("output/MNIST_labels_70k.npy")

print(dim(X))

CORES <- 30
REPLICATIONS <- 2

## Implement cluster Accuracy equivalent to Python implementation
cluster_accuracy <- function(y_true, y_pred) {
  cm <- as.matrix(table(y_true, y_pred))
  if (nrow(cm) > ncol(cm)) {
    cm <- t(cm)
  }
  assignment <- solve_LSAP(cm, maximum = TRUE)
  matched_count <- sum(cm[cbind(seq_along(assignment), assignment)])
  return(matched_count / length(y_true))
}


## PHM-Ensemble 
cat("Starting Globalized PHM-Ensemble!\n")
SUBSAMPLE_SIZE <- 7000
for (idx in 1:REPLICATIONS) {
  elapsed_all <- system.time({
    pars3 <- constructPmcParamsSubAggMclust(X,
                                            50, SUBSAMPLE_SIZE,
                                            G=1:50,
                                            numCores=CORES,
                                            modelNames=c("VVV"),
                                            prefix=paste0("mclustSubaggALL", num_pcs, "_r", idx),
                                            saveDir=ensemble_directory,
                                            verbose=F)    
  })
  all_time_fname <- file.path(ensemble_directory, paste0("mclustSubaggALL", num_pcs, "_r", idx, "_timing.npy"))
  if (!file.exists(all_time_fname)) {
    np$save(all_time_fname,
            unname(elapsed_all[3]))
  }

  cat("\tFinished Replicate", idx, ":",  elapsed_all[3], "seconds\n")
}


## Single-Replicate PHM based on hierarchical initialization
cat("Starting Single-Replicate PHM!\n")
temp <- mclapply(1:REPLICATIONS, function(idx) {
  fname <- file.path(single_directory, paste0("mclustALL", num_pcs, "_G50", "_r", idx, ".RDS"))
  if (!file.exists(fname)) {
    elapsed <- system.time({
      set.seed(idx)
      mcl4 <- Mclust(X, G=1:50, modelNames=c("VVV"), verbose=F)
      pars4 <- constructPmcParamsMclust(mcl4)
    })
    saveRDS(pars4, file=fname)
    np$save(file.path(single_directory, paste0("mclustALL", num_pcs, "_G50", "_r", idx, "_timing.npy")),
            unname(elapsed[3]))
  }
  return(idx)
}, mc.cores=CORES)
cat("Finished Single-Replicate PHM!\n")


## Mclust (GMM w/ hierarchical initialization) with G=10 fixed
cat("Starting mclust w/ G=10!\n")
temp <- mclapply(1:REPLICATIONS, function(idx) {
  fname_params <- file.path(mclust_directory, paste0("mclustG10_params_r", idx, ".RDS"))
  fname_pred <- file.path(mclust_directory, paste0("mclustG10_pred_r", idx, ".npy"))
  fname_time <- file.path(mclust_directory, paste0("mclustG10_timing_r", idx, ".npy"))
  
  if (file.exists(fname_pred)) return(idx)
  elapsed <- system.time({
    set.seed(idx)
    mclSingle <- Mclust(X, G=10)
  })
  np$save(fname_pred, mclSingle$classification)
  np$save(fname_time, unname(elapsed[3]))
  paramsSingle <- constructPmcParamsMclust(mclSingle)
  saveRDS(paramsSingle, file=fname_params)
  return(idx)
}, mc.cores=CORES)
cat("Finished mclust w/ G=10!\n")


################################
### Oracle Density Estimator ###
################################
cat("Starting Oracle Density!\n")
if (file.exists(gt_file)) {
    load(gt_file)
} else {
    set.seed(1)
    pars_gt = list()
    for (lab in 0:9) {
      print(lab)
      lab_idx = which(labels == lab)
      sub_X = X[lab_idx, ]
    
      class_weight = nrow(sub_X) / nrow(X)
      sub_pars = constructPmcParamsSubAggMclust(sub_X, 25, 3000,
                                                G=1:20, numCores=CORES)
    
      new_pars = list(
        prob=numeric(length(sub_pars)),
        mean=matrix(nrow=num_pcs, ncol=length(sub_pars)),
        var=array(dim=c(num_pcs, num_pcs, length(sub_pars))),
        class=lab
      )
      for (idx in 1:length(sub_pars)) {
        new_pars$prob[idx] <- sub_pars[[idx]]$prob * class_weight
        new_pars$mean[, idx] <- sub_pars[[idx]]$mean[, 1]
        new_pars$var[, , idx] <- sub_pars[[idx]]$var[, , 1]
      }
    
      pars_gt[[lab+1]] <- new_pars
    }
    
    deltapmc <- computeMonteCarloDeltaPmcMatrix(pars_gt, 
                                                mcSamples=1e6, 
                                                batchSize=1e4,
                                                numCores=CORES,
                                                verbose=T)
    
    save(pars_gt, deltapmc, 
         file=gt_file)
}

post_mat <- computePosteriorProbMatrix(pars_gt, X)
oracle_labels <- apply(post_mat, 1, which.max)
cat("Completed Oracle Density!\n")


###############################
### Hierarchical Clustering ###
###############################
cat("Started Hierarchical Clustering!\n")
elapsed_hcl <- system.time({
  hcl <- fastcluster::hclust(
    dist(X), method="ward.D2"
  )  
})

hclust_labels <- cutree(hcl, 10)
cat("Completed Hierarchical Clustering!\n")

print(knitr::kable(data.frame(
    Method=c("Oracle", "Hier. Clust"),
    ARI=c(
        ARI(factor(labels), factor(oracle_labels)),
        ARI(factor(labels), factor(hclust_labels))
    ),
    NMI=c(
        NMI(factor(labels), factor(oracle_labels)),
        NMI(factor(labels), factor(hclust_labels))
    ),
    Acc=c(
        cluster_accuracy(labels, oracle_labels),
        cluster_accuracy(labels, hclust_labels)
    )
), digits=3))