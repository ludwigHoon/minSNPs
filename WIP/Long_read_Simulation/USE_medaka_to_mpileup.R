read_vcf <- function(file) {
    read_data <- readLines(file)
    header_line <- max(which(startsWith(read_data, "#")))
    header <- gsub("#", "", read_data[header_line])
    data <- paste(c(header, read_data[(header_line + 1):length(read_data)]), collapse = "\n")
    result <- read.table(text = data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(result)
}

filter_to_snp <- function(vcf_data) {
    is_snp <- nchar(vcf_data$REF) == 1 & nchar(vcf_data$ALT) == 1
    return(vcf_data[is_snp,])
}

get_snp_pos_from_vcfs <- function(vcfs, barcodes) {
    outputs <- list()
    if (is.null(barcodes)) {
        barcodes <- seq_along(vcfs)
    } else {
        stopifnot(length(vcfs) == length(barcodes))
    }
    for (v in seq_along(vcfs)) {
        vcf_data <- read_vcf(vcfs[v])
        snp_data <- filter_to_snp(vcf_data)
        snp_data$barcode <- barcodes[v]
        outputs[[as.character(v)]] <- snp_data
    }
    all_snp_data <- rbindlist(outputs)
    common_snps <- Reduce(intersect, sapply(outputs, `[[`, "POS"))
    return(list(
        all_snp_data = all_snp_data,
        common_snps = common_snps
    ))
}

all_vcf <- list.files(pattern = "medaka.vcf", path="sup", recursive = TRUE, full.names = TRUE)
barcodes <-  sapply(strsplit(all_vcf, split = "_"), `[`, 2)
all_snp_data <- get_snp_pos_from_vcfs(all_vcf, barcodes)

CHROM <- unique(all_snp_data$all_snp_data$CHROM)
POS <- unique(all_snp_data$all_snp_data$POS)

stopifnot(length(CHROM) == 1)

result <- data.frame(chrom = rep(CHROM, length(POS)), pos = sort(POS, decreasing = FALSE))
write.table(result, file = "SNPS_POS_PILEUP.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

ref_fasta <- "/home/lhoon/MINSNPS_P2/Mu50.fasta"
snps_pos_pileup <- "SNPS_POS_PILEUP.txt"
medaka_calls_to_ref <- list.files(pattern = "calls_to_ref.bam", path="sup", recursive = TRUE, full.names = TRUE)


writeLines("#!/bin/bash

#SBATCH --job-name=medaka_haploid_variant   ## name that will show up in the queue
#SBATCH --ntasks=24  ## number of tasks (analyses) to run
#SBATCH --cpus-per-task=2  ## the number of threads allocated to each task
#SBATCH --time=1-00:00:00  ## time for analysis (day-hour:min:sec)

file=`ls --directory sup/*_medaka`;

for f in $file;
do
  srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK bash -c \"conda activate medaka; bcftools mpileup --threads 2 --skip-indels --count-orphans --no-BAQ --min-BQ 0 -R SNPS_POS_PILEUP.txt -Ov --annotate DP,AD -f /home/lhoon/MINSNPS_P2/Mu50.fasta -o $f/all_snp.vcf $f/calls_to_ref.bam\" &
done

wait", "bcftools_mpileup.sh")
system("sbatch bcftools_mpileup.sh")