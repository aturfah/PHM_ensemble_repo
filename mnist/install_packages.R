options(repos = c("https://cloud.r-project.org/"))

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

if (!require("PHM", quietly = TRUE)) {
    remotes::install_github("aturfah/PHM")
}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("rhdf5", quietly = TRUE))
    BiocManager::install("rhdf5")

if (!require("aricode", quietly = TRUE))
    install.packages("aricode")

if (!require("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
}

if (!require("parallel", quietly = TRUE))
    install.packages("parallel")

if (!require("clue", quietly = TRUE))
    install.packages("clue")

if (!require("fastcluster", quietly = TRUE))
    install.packages("fastcluster")

