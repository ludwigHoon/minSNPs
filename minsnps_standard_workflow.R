# MinSNPs package standard workflow
## Instructions:
# 1. Edit the "EDIT THIS" section with relevant info
# 2. Run Rscript minsnps_standard_workflow.R

#################################################################### EDIT THIS ####################################################################

###################### Input Options ######################
input_alignment <- "./Orthologous_SNP_Matrix.fasta"
metadata <- "./metadata.csv" # Must contain at least 2 columns - isolate, target
# Only isolate in both metadata and input alignment will be used (must have all the sequences names in the input alignment)

###################### Multiprocessing Options ######################
# Default to "default" <the number of cores detected - 2>, set to 1 or more cores to use 
nworkers <- 2

###################### read_fasta Options ######################
## read alignment option, either TRUE or FALSE
force_to_upper <- TRUE

###################### process_allele Options ######################
# Check alignment option
check_length <- TRUE # if TRUE, will check the length of each sequence alignment and remove those longer/shorter than most other isolates
check_bases <- TRUE # if TRUE, check each bases, with the following settings
### Will only matter if check_bases is TRUE
dash_ignore <- TRUE # If TRUE, positions containing "-" will be ignored, otherwise treated as a different type
accepted_char <- c("A", "C", "T", "G") # use upper character, if read_fasta force matrix to upper case is set to true
ignore_case <- TRUE # If TRUE, both alignment and accepted_char will be force-converted to upper case
remove_invariant <- TRUE # If TRUE, remove positions where there is only 1 allele
biallelic_only <- FALSE # If TRUE, remove positions with more than 2 alleles


###################### find_optimised_snps Options ######################
# Available modes are: 
# - Discriminate some from all:
#     - "percent" (Constrained to 100% sensitivity)
#     - "mcc_single" (Binary-class MCC) # Version 0.2.0 onwards
# - Discriminate all from all:
#     - "simpson" (Based on Simpson's index of diversity)
#     - "mcc_multi" (Multi-class MCC) # Version 0.2.0 onwards
#     - "simpson_by_group" (Averaging of Simpson's index of diversity and Multi-class MCC) # Version 0.2.0 onwards
#
# 
analysis_mode <- "percent" # Available modes: "simpson", "percent"; version >= 0.2.0 also supports mcc, mcc_multi, simpson_by_group
number_of_result <- 5 # number of sets
max_depth <- 5 # maximum number of SNPs per set
goi_target <- "SEA" # or suitable tag in the target column of metadata
accept_multiallelic <- TRUE # For `Discriminate some from all` mode, whether the positions can contain variants position within group of interest
included_positions = c()
excluded_positions = c()
search_from = NULL # if not NULL, all positions not listed here are excluded, otherwise use a vector `c()`
iterate_included = FALSE # If TRUE, output the score of the included_positions at each level
completely_unique = FALSE # If TRUE, all the SNPs returned in each set will be completely different, default to only requiring the first SNPs to be different (This will be slightly faster when requiring number_of_result > 1)
priority <- NULL # Used for mcc_multi and simpson_by_group,
# prioritise tagging allele profile to class in descending order
# If NULL, will be automatically computed based on target column,
# target with the least number of sample will be prioritised first,
# tie-breaking with target value alphabetically.
# To supply, all targets value must be in the vector, e.g.,
# c("Target_val_1", "Target_val_2", ....)
# where Target_val_1 is prioritised over Target_val_2

###################### output_result Options ######################le 
output_filename <- "minSNPs_output.tsv" # Output file name, output is a tab-delimited file

#################################################################### END ####################################################################
############################################## Only touch below if you know what you are doing ##############################################
# Libraries
if (!require("minSNPs", quietly = TRUE)){
  if (!require("BiocParallel", quietly = TRUE)){
    if (!require("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("BiocParallel")
  }
  install.packages("minSNPs")
}

library(minSNPs)
library(BiocParallel)

timestamp()
print("Analysis started")

# Good options
options(stringsAsFactors = FALSE)

meta <- read.csv(metadata)
stopifnot(c("isolate", "target") %in% colnames(meta))
stopifnot( (is.numeric(nworkers) || nworkers == "default") )

if (nworkers == 1){
    bp <- SerialParam()
} else {
    if (nworkers == "default"){
        nworkers <- multicoreWorkers()
    }
    bp <- MulticoreParam(workers = as.numeric(nworkers))
}

print(paste0("Using ", nworkers, " cores for processing"))

# Read SNP Matrix
snp_matrix <- read_fasta(input_alignment, force_to_upper = force_to_upper)

isolate_in_matrix <- names(snp_matrix)
isolate_in_meta <- meta$isolate
common_isolate <- intersect(isolate_in_matrix, isolate_in_meta)

if (length(common_isolate) <= 1)
    stop("Only 1 or less common isolates found in matrix")

if (length(common_isolate) != length(isolate_in_matrix) || length(common_isolate) != length(isolate_in_meta)){
    warning("Not complete match of isolates in metadata and input alignment, only subset of these are used")
}

meta <- meta[meta$isolate %in% common_isolate,]
snp_matrix <- snp_matrix[common_isolate]

# Preprocess SNP Matrix
p_snp_matrix <- process_allele(seqc = snp_matrix, bp = bp, check_length=check_length,
                               check_bases=check_bases, dash_ignore=dash_ignore,
                               accepted_char=accepted_char, ignore_case=ignore_case,
                               remove_invariant=remove_invariant, biallelic_only=biallelic_only)

goi <- as.character(meta[meta$target == goi_target,]$isolate)

if (any(!unique(analysis_mode) %in% names(get_metric_fun()))){
    stop("Some analysis modes are not supported; NOTE: minSNPs prior to v0.2.0 supports only percent and simpson modes")
}

### version under 0.2.0 does not support mode other than "simpson" and "percent"
if (unlist(packageVersion("minSNPs"))[2] >= 2){
    if (is.null(priority))
        priority <- generate_prioritisation(meta)
}


result <- find_optimised_snps(seqc = p_snp_matrix, metric = analysis_mode, goi = goi,
    accept_multiallelic = accept_multiallelic, number_of_result = number_of_result, max_depth = max_depth,
    included_positions = included_positions, excluded_positions = excluded_positions, search_from = search_from,
    iterate_included = iterate_included, completely_unique = completely_unique, bp = bp, output_progress = TRUE,
    target = "target", meta = meta, priority = priority)
output_result(result, view = "tsv", file_name = output_filename)

timestamp()
print("Analysis completed")
