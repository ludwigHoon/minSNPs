### HELPERS
STEP_1 <- T

process_result_file <- function(result_filepath){
    f <- file(result_filepath, "r")
    preparse <- readLines(f)
    result <- list()
    cur <- 1
    is_result <- F
    for (i in seq_len(length(preparse))){
        if (length(grep("^Result", preparse[[i]])) >= 1){
            is_result <- T
        }
        if (is_result){
            if (preparse[[i+1]] == ""){
                result[[cur]] <- as.numeric(strsplit(gsub("\"", "", strsplit(preparse[[i]], "\t")[[1]][1]), ", ")[[1]])
                cur <- cur + 1
                is_result <- F
            }
        }
    }
    close(f)
    final_selected <- (result)
    return(final_selected)
}

generate_kmers <- function(final_string, k){
    if (is.character(final_string) && length(final_string) == 1){
        final_string <- strsplit(final_string, split = "")[[1]]
    }
    kmer <- list()
    for (i in seq_len(length(final_string) - k + 1)){
        kmer[[i]] <- paste(final_string[i:(i+k-1)], collapse = "")
    }
    return(unlist(kmer))
}

process_variant_file <- function(variant_sites){
    result <- c()
    read_cc <- c()
    con <- file(variant_sites, "r")
    while(TRUE){
        data <- readLines(con, n = 3)
        if (length(data) == 0){
            break
        }
            positions <- strsplit(strsplit(data, split = "\n")[[2]], split = ", ")
            read_cc <- c(read_cc, strsplit(data, split = "\n")[[1]])
            result <- c(result, unlist(positions))
    }
    close(con)
    print(read_cc)
    result <- sort(as.numeric(unique(result)))
    return(result)
}

reverse_complement <- function(seq){
    seq <- strsplit(seq, split = "")[[1]]
    result <- lapply(seq, function(x){
        if (toupper(x) == "A") {return("T")}
        if (toupper(x) == "C") {return("G")}
        if (toupper(x) == "T") {return("A")}
        if (toupper(x) == "G") {return("C")}
        return(x)
        })
    return(rev(paste(result, collapse="")))
}



### STEP 1
identify_overlaps <- function(selected_fasta, position_reference, prev, after, extend_length = T){
    snp_table <- data.frame(snp_id = c(), fasta_position = c(), genome_position = c(), string_start = c(), string_end = c(),
        stringsAsFactors = F)
    overlap_table <- data.frame(snp_id = c(), overlaps_genome_pos = c(), overlaps_fasta_pos = c(),
        overlaps_string_pos = c(), stringsAsFactors = F)

    selected_fasta <- c(selected_fasta)
    for (snp in selected_fasta){
        genome_pos <- position_reference[position_reference$fasta_position == snp, "genome_position"]
        string_start <- genome_pos - prev
        string_end <- genome_pos + after
        
        N_test <- 0
        snp_string_pos <- prev + 1
        # Check overlaps before
        overlaps_before <- position_reference[position_reference$genome_position >= (genome_pos - (N_test+prev)) &
            position_reference$genome_position < genome_pos, c("fasta_position", "genome_position")]
        if(extend_length){
            while((N_test + prev) - nrow(overlaps_before) < prev){
                N_test <- N_test + (nrow(overlaps_before) - N_test)
                overlaps_before <- position_reference[position_reference$genome_position >= (genome_pos - (N_test+prev)) &
                    position_reference$genome_position < genome_pos, c("fasta_position", "genome_position")]
            }
            string_start <- string_start - N_test
            snp_string_pos <- snp_string_pos + N_test
            overlaps_before$s_pos <- overlaps_before$genome_position - string_start
        }
        
        N_test <- 0
        # Check overlaps after
        overlaps_after <- position_reference[position_reference$genome_position > genome_pos &
            position_reference$genome_position <= (genome_pos + (N_test+after)), c("fasta_position", "genome_position")]
        if(extend_length){
            while((N_test + after) - nrow(overlaps_after) < after){
                N_test <- N_test + (nrow(overlaps_after) - N_test)
                overlaps_after <- position_reference[position_reference$genome_position > genome_pos &
                    position_reference$genome_position <= (genome_pos + (N_test+after)), c("fasta_position", "genome_position")]
            }
            string_end <- string_end + N_test
            overlaps_after$s_pos <- snp_string_pos + (overlaps_after$genome_position - genome_pos)
        }
        
        temp_overlap_table <- rbind(overlaps_before, overlaps_after)
        if (ncol(temp_overlap_table) >= 3){
            colnames(temp_overlap_table) <- c("overlaps_fasta_pos", "overlaps_genome_pos", "overlaps_string_pos")
        }
        temp_overlap_table$snp_id <- rep(snp, nrow(temp_overlap_table))

        snp_table <- rbind(snp_table, data.frame(
            snp_id = snp, fasta_position = as.character(paste(snp, collapse = ", ")),
                genome_position = genome_pos, string_start = string_start,
                string_end = string_end, snp_string_pos = as.character(paste(snp_string_pos, collapse = ", ")),
                stringsAsFactors = F))
        overlap_table <- rbind(overlap_table, temp_overlap_table)
    }
    return(list(snp_table = snp_table, overlap_table = overlap_table))
}

#data.frame(snp_id = c(), genome_position = c(), fasta_position = c(), string_start = c(),
#    string_end = c(), snp_string_pos = c(), stringsAsFactor = F)

#data.frame(snp_id = c(), overlaps_genome_pos = c(), overlaps_fasta_pos = c(),
#    overlaps_string_pos = c(), stringsAsFactor = F)
match_count <- function(target, search_from){
    matches <- gregexpr(paste(target, collapse=""), search_from)
    found <- 0
    if (length(matches[[1]]) == 1){
        if (matches[[1]][1] != -1){
            found <- 1
        }
    } else{
        found <- length(matches[[1]])
    }
    return(found)
}

### search string for gene
generate_search_string_gene <- function(genes, ref_seq, k, id_prefix="s"){
    
    string_table <- data.frame(search_string = c(), string_id = c(), stringsAsFactors = F)

    mapping_table <- data.frame(gene = c(), f_string_id = c(), r_string_id = c(), 
        n_match_genome = c(), n_match_genome_rev = c(), stringsAsFactors = F)

    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }
    
    search_strings <- BiocParallel::bplapply(genes, function(gene, generate_kmers, k){
        kmers <- generate_kmers(final_string = gene, k = k)
        return(kmers)
    }, k=k, generate_kmers=generate_kmers)

    names(search_strings) <- names(genes)

    u_string <- unique(unlist(search_strings))
    rc_u_string <- unname(sapply(u_string, reverse_complement))

    string_table <- data.frame(search_string = c(u_string, rc_u_string),
        string_id = c(paste(id_prefix, seq_along(u_string), sep = "_"),
            paste(id_prefix, (seq_along(rc_u_string) + length(u_string)), sep = "_")),
        stringsAsFactors = F)

####HERE
    string_table$ 

    temp_mapping_list <- BiocParallel::bplapply(names(search_strings),
        function(i, search_strings, string_table, id_prefix){
            search_string <- search_strings[[i]]
        temp_mapping_table <- data.frame(gene = c(), f_string_id = c(), r_string_id = c(), 
            stringsAsFactors = F)
        
        string_id <- match(search_string, string_table$search_string)
        f_string_id <- paste(id_prefix, string_id, sep = "_")
        r_string_id <- paste(id_prefix, (string_id + nrow(string_table)/2), sep = "_")
        return(data.frame(gene = rep(i, length(string_id)),
            f_string_id = f_string_id, r_string_id = r_string_id, stringsAsFactors = F))
    }, string_table=string_table, search_strings=search_strings, id_prefix=id_prefix)

    mapping_table <- do.call(rbind, temp_mapping_list)
    
    data.frame(gene = names(genes), f_string_id = c(), r_string_id = c(), 
        n_match_genome = c(), n_match_genome_rev = c(), stringsAsFactors = F)


    temp_string_table$n_match_genome <- unlist(lapply(temp_string_table$search_string, match_count,
        search_from = paste(ref_seq, collapse = "")))
    temp_string_table$n_match_genome_rev <- unlist(lapply(temp_string_table$search_string, match_count,
        search_from = reverse_complement(paste(ref_seq, collapse = ""))))



    search_string <- c()
    for (gene in genes){
        gene_seq <- ref_seq[gene]
        gene_seq <- strsplit(gene_seq, split = "")[[1]]
        gene_seq <- paste(gene_seq, collapse = "")
        gene_seq <- paste(gene_seq, reverse_complement(gene_seq), sep = "|")
        search_string <- c(search_string, gene_seq)
    }
    search_string <- paste(search_string, collapse = "|")
    search_string <- paste0("(?<=.{", k, "})", search_string, "(?=.{", k, "})")
    
    return(search_string)
}

#data.frame(search_string = c(), snp_id = c(), genes = c(), snp_string = c(),
#    n_match_genome = c(), n_match_genome_rev = c(), stringsAsFactor = F)



### STEP 2
generate_search_string <- function(snp_table, overlap_table, orth_matrix, ref_seq, include_neighbour = F,
    bpparam = BiocParallel::MulticoreParam()){

    string_table <- data.frame(search_string = c(), snp_id = c(), snp_string = c(),
        n_match_genome = c(), n_match_genome_rev = c(), stringsAsFactors = F)
    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }
    
    #pb = txtProgressBar(min = 0, max = nrow(snp_table), initial = 0, style = 3) 
    l_string_table <- BiocParallel::bplapply(seq_len(nrow(snp_table)), function(i){
        #setTxtProgressBar(pb,i)
        snp_id <- snp_table[i, "snp_id"]
        #tic(snp_id)
        ref_string <- ref_seq[c(snp_table[i, "string_start"]:snp_table[i, "string_end"])]
        fasta_positions <- as.numeric(strsplit(snp_table[i, "fasta_position"], split = ", ")[[1]])
        overlaps <- overlap_table[overlap_table$snp_id == snp_id,]

        if (include_neighbour){
            ##### TO-DO
            next
        } else {
            variants <- unlist(unique(minSNPs:::generate_pattern(orth_matrix, fasta_positions)))
            f_string <- lapply(variants, function(var){
                ref_string[as.numeric(strsplit(snp_table[i, "snp_string_pos"], split = ", ")[[1]])] <- var
                return(paste(ref_string, collapse = ""))
            })
            f_string <- lapply(f_string, function(string){
                string <- strsplit(string, split = "")[[1]]
                string[overlaps$overlaps_string_pos] <- "."
                return(paste(string, collapse = ""))
            })
        }
        rc_string <- lapply(f_string, reverse_complement)
        temp_string_table <- data.frame(search_string = c(unlist(f_string), unlist(rc_string)),
            snp_id = rep(snp_id, length(f_string)*2), snp_string = c(variants, variants), stringsAsFactors = F)
        temp_string_table$n_match_genome <- unlist(lapply(temp_string_table$search_string, match_count,
            search_from = paste(ref_seq, collapse = "")))
        temp_string_table$n_match_genome_rev <- unlist(lapply(temp_string_table$search_string, match_count,
            search_from = reverse_complement(paste(ref_seq, collapse = ""))))
        temp_string_table$fasta_position <- rep(snp_table[i, "fasta_position"], nrow(temp_string_table))
        #toc()
        return(temp_string_table)
    }, BPPARAM = bpparam)
    string_table <- do.call(rbind, l_string_table)
    #close(pb)
    return(string_table)
}

#data.frame(search_string = c(), snp_id = c(), fasta_position = c(), snp_string = c(),
#    n_match_genome = c(), n_match_genome_rev = c(), stringsAsFactor = F)


### STEP 3
summarise_result <- function(snp_table, overlap_table, search_string_table, excluded_position = NULL){
    new_snp_table <- snp_table
    snps <- unlist(unique(snp_table$snp_id))
    new_snp_table$multiple_match <- unlist(lapply(snps, function(id){
        subset <- snp_table[snp_table$snp_id == id, ]
        return((sum(subset$n_match_genome) > 1) | (sum(subset$n_match_genome_rev) > 1))
    }))
    new_snp_table$variant_count <- unlist(lapply(snps, function(id){
        subset <- search_string_table[search_string_table$snp_id == id,]
        return(nrow(subset)/2)
    }))
    new_snp_table$contain_excluded <- unlist(lapply(snps, function(id){
        subset <- overlap_table[overlap_table$snp_id == id, "overlaps_fasta_pos"]
        return(any(subset %in% excluded_position))
    }))
    return(new_snp_table)
}



#### RUNS
library(minSNPs) 
library(BiocParallel) 

ref_matrix <- read.csv("ref_OUTF1.csv")
ort_matrix <- read_fasta("./balanced_mat_without_SRR798725_single_final.fasta")
excluded <- process_variant_file("variant_t10.csv")
result_filepath <- "evenhd_with_exc200x8_400.csv"
final_selected <- unique(unlist(process_result_file(result_filepath)))
ref_genome <- read_fasta("./Mu50.fasta")

###
# final_selected <- seq_len(length(ref_matrix[[1]]))
# final_selected <- final_selected[!final_selected %in% excluded]
##
prev <- 7
after <- 7

if (STEP_1) {
RES_STEP_1 <- bplapply(final_selected, identify_overlaps,
    position_reference = ref_matrix, prev = prev, after = after)

snp_table <- do.call(rbind, lapply(RES_STEP_1, function(x){
    return(x$snp_table)}))
overlap_table <- do.call(rbind, lapply(RES_STEP_1, function(x){
    return(x$overlap_table)}))

write.csv(snp_table, "snp_table.csv", row.names = F)
write.csv(overlap_table, "overlap_table.csv", row.names = F)

RES_STEP_1[["snp_table"]] <- snp_table
RES_STEP_1[["overlap_table"]] <- overlap_table
}

search_string_table <- generate_search_string(RES_STEP_1[["snp_table"]],
                        RES_STEP_1[["overlap_table"]], ort_matrix, ref_genome, bpparam = MulticoreParam(workers = 6, progressbar = T))
write.csv(search_string_table, "search_string_table.csv", row.names = F)

summary <- summarise_result(RES_STEP_1[["snp_table"]], RES_STEP_1[["overlap_table"]],
                        search_string_table)
write.csv(summary, "summary.csv", row.names = F)

#search_string_table <- generate_search_string(test_snp_table, test_overlap_table, ort_matrix, ref_genome, bpparam = SerialParam(progressbar = T))