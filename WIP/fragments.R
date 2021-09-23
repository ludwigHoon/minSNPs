#' \code{read_ref_meta}
#'
#' @description
#' \code{read_ref_meta} is used to read the reference files (genome and meta).
#' @param ref_meta path to the meta file
#' @param ref_genome either path to the refence genome in fasta format
#' or the read sequence (e.g. from \code{read_fasta}).
#' @param read_genome_from_file whether to read the reference genome from file
#' @param read_meta_from_file whether to read the meta from file
#' @return Will return a list consist of the genome and a dataframe of the meta
read_ref_meta <- function(ref_meta, ref_genome, read_meta_from_file = TRUE,
    read_genome_from_file = TRUE) {
    if (read_meta_from_file) {
        meta <- read.csv(ref_meta)
    } else{
        meta <- ref_meta
    }
    if (!"fasta_position" %in% names(meta) ||
        !"genome_position" %in% names(meta)) {
            stop("Both fasta_position and genome position may be present")
    }
    if (read_genome_from_file) {
        genome <- minSNPs::read_fasta(ref_genome)
        if (length(genome) > 1) {
            cat("More than 1 genome read,",
            "using the 1st one, ignoring the rest\n")
            genome <- genome[[1]]
        }
        if (length(genome) < 1) {
            stop("Can't read any genome\n")
        }
    } else {
        genome <- ref_genome
    }
    result <- list(meta = meta, genome = genome)
    class(result) <- "ref_meta_genome"
    return(result)
}

get_ref_genome_position <- function(fasta_positions = c(), meta_df) {
    positions <- c(fasta_positions)
    if (class(meta_df) == "ref_meta_genome") {
        meta_df <- meta_df$meta
    }
    ref_meta <- meta_df
    result <- ref_meta[match(positions, ref_meta$fasta_position),
        "genome_position"]
    names(result) <- positions
    not_found_positions <- names(result[is.na(result)])
    if (length(not_found_positions) > 0) {
        cat("These Fasta Position is not mapped to reference genome",
        "and is ignored", paste(not_found_positions, collapse = ", "), "\n")
    }
    return(result[!is.na(result)])
}

# regions <-
#   [[region1]]
#   [[region1]][[start]]
#   [[region1]][[end]]
#   [[region2]]
#   [[region2]][[start]]
#   [[region2]][[end]]
generate_fragments <- function(regions = c(), sequence) {
    regions <- c(regions)
    fragments <- bplapply(regions, function(reg, sequence) {
        start <- reg[["start"]]
        end <- reg[["end"]]
        reg[["sequence"]] <- sequence[[1]][start:end]
        return(reg)
    }, sequence = sequence)
    return(fragments)
}

collapse_regions <- function(regions) {
    # First order the region according to starting position
    ordered_regs <- regions[order(sapply(regions, `[[`, "start"))]
    new_regions <- list()
    frag_i <- 1
    skip <- c()
    # Go through each regions
    for (i in seq_len(length(ordered_regs))) {
        # Do not have to go through this region,
        # since we've added it to skip list
        if (i %in% skip) {
            next
        }
        # This is the last region, but not skipped,
        # So we add it to the end
        if (i == length(ordered_regs)) {
            reg <- list()
            reg[["start"]] <- ordered_regs[[i]][["start"]]
            reg[["end"]] <- ordered_regs[[i]][["end"]]
            new_regions[[paste("region-", frag_i, sep = "")]] <- reg
            break
        }
        # Otherwise, we create the starting position of the region
        reg <- list()
        reg[["start"]] <- ordered_regs[[i]][["start"]]

        # If the ending position of this region is >=
        # the starting of the next one,
        if (ordered_regs[[i]][["end"]] >= ordered_regs[[i + 1]][["start"]]) {
            # We first take the ending position of the next one
            reg[["end"]] <- ordered_regs[[i + 1]][["end"]]
            # And iteratively check if next ending position is >=
            # the current one
            for (j in seq.int((i + 1), length(ordered_regs))) {
                # This is the last one, there will be no next one to check
                if (j == length(ordered_regs)) {
                    break
                }
                # If the next starting < current ending
                if (reg[["end"]] >= ordered_regs[[j + 1]][["start"]]) {
                    reg[["end"]] <- ordered_regs[[j + 1]][["end"]]
                    skip <- c(skip, j + 1)
                } else {
                    # No need to continue checking, since broken
                    break
                }
            }
            skip <- c(skip, i + 1)
        } else {
            # Reassigning the ending position, no modification needed
            reg[["end"]] <- ordered_regs[[i]][["end"]]
        }
        new_regions[[paste("region-", frag_i, sep = "")]] <- reg
        frag_i <- frag_i + 1
    }
    return(new_regions)
}

## TODO
identify_mutations <- function(positions = c(), reference, alignment) {
    positions <- c(positions)
    result <- list()
    return(result)
}

## TODO
generate_meta <- function(regions = c(), sequence) {
    regions <- c(regions)
    meta <- bplapply(regions, function(reg, sequence) {
        start <- reg[["start"]]
        end <- reg[["end"]]
        return(sequence[[1]][start:end])
        }, sequence = sequence
    )
    return(meta)
}

get_snp_in_region <- function(meta_df, regions = c()) {
    regions <- c(regions)
    meta <- meta_df$meta
    result <- bplapply(regions, function(reg, meta) {
        all_snps <- meta[meta["genome_position"] <= reg[["end"]] &
            meta["genome_position"] >= reg[["start"]],
                c("genome_position", "fasta_position")]
        reg[["SNPs_ref_genome"]] <- all_snps[["genome_position"]]
        reg[["SNPs_fasta"]] <- all_snps[["fasta_position"]]
        reg[["SNPs_frag"]] <- reg[["SNPs_ref_genome"]] - reg[["start"]] + 1
        return(reg)
    }, meta = meta)
    return(result)
}


get_fragments_isolates <- function(meta_df, fasta_positions = c(),
    before_snp_length = 150, after_snp_length = 150) {

    genome_length <- length(meta_df$genome[[1]])
    fasta_positions <- c(fasta_positions)

    # These are the position of interested SNP in reference genome
    genome_positions <- get_ref_genome_position(fasta_positions, meta_df)

    # Check the start and end, force to 1 or length of sequence if necessary
    regions <- bplapply(genome_positions, function(pos) {
        return(list(start =
            ifelse((pos - before_snp_length) < 1, 1, (pos - before_snp_length)),
        end = ifelse((pos + after_snp_length) > genome_length,
            genome_length, (pos + after_snp_length))
        ))
    })

    # Then collapse the overlapping regions
    regions <- collapse_regions(regions)

    # Position of other SNPs in reference genome within the regions
    oth_genome_positions <- get_snp_in_region(meta_df, regions)

    # Generate the fragments
    result <- generate_fragments(oth_genome_positions, meta_df$genome)
    return(result)
}


## Position in Fragment = SNPs_ref_genome - start + 1 <<<----- FOR META GENERATION


# list of fragments
# fragments <-
# [[fragment 1]]
# [[fragment 1]][[start]] = 10  (reference genome position)
# [[fragment 1]][[end]] = 30 (reference genome position)
# [[fragment 1]][[sequence]] = "A C TG"
# [[fragment 1]][[SNPs - fasta]] = c(12,....)
# [[fragment 1]][[SNPs - ref]] = c(12,....)
# [[fragment 1]][[wild_type]] = "ACTG...." * >= 10 to <= 30
# [[fragment 1]][[mutated]] = c("ACTG....", "ACTG....")


# [[fragment 2]]
# [[fragment 2]][[genome_start]] = 10
# [[fragment 2]][[genome_end]] = 30
# [[fragment 2]][[SNPs]] = 12
# [[fragment 2]][[wild_type]] = "ACTG...." * >= 10 to <= 30
# [[fragment 2]][[mutated]] = c("ACTG....", "ACTG....")


# [[fragment 3]]
# [[fragment 3]][[genome_start]] = 10
# [[fragment 3]][[genome_end]] = 30
# [[fragment 3]][[SNPs]] = 12
# [[fragment 3]][[wild_type]] = "ACTG...." * >= 10 to <= 30
# [[fragment 3]][[mutated]] = c("ACTG....", "ACTG....")