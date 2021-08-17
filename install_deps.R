if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocParallel")
BiocManager::install("Biostrings")