### Generate search string for gene

generate_search_string_gene <- function(gene_seq, ref_seq, k, id_prefix = "gene"){
    string_table <- data.frame(search_string = c(), string_id = c(), stringsAsFactors = FALSE)

    if (class(ref_seq) == "list") {
        ref_seq <- ref_seq[[1]]
    }

    search_strings <- BiocParallel::bplapply(gene_seq, function(gene, generate_kmers, k){
        kmers <- generate_kmers(final_string = gene, k = k)
        return(kmers)
    }, k=k, generate_kmers=generate_kmers)

    u_string <- unique(unlist(search_strings))
    rc_u_string <- unname(sapply(u_string, reverse_complement))

    string_table <- data.frame(search_string = c(u_string, rc_u_string),
        string_id = c(paste(id_prefix, seq_along(u_string), "1", sep = "_"),
            paste(id_prefix, seq_along(rc_u_string), "2", sep = "_")),
        stringsAsFactors = F)

    string_table$n_match_genome <- unlist(bplapply(string_table$search_string, match_count,
        search_from = paste(ref_seq, collapse = "")))
    string_table$n_match_genome_rev <- unlist(bplapply(string_table$search_string, match_count,
        search_from = reverse_complement(paste(ref_seq, collapse = ""))))

    return(string_table)
}

collapse_strings_table <- function(string_table_1, string_table_2,
    combine_to = data.frame(search_strings = c(), string_id = c(),
        stringsAsFactors = FALSE)){

    all_strings <- c(string_table_1$search_string, string_table_2$search_string,
        combine_to$search_strings)
    all_strings <- unique(unlist(all_strings))

    collapsed_strings_table <- data.frame(search_string = all_strings)

    strings_id <- lapply(collapsed_strings_table$search_string,
        function(x){
            collapsed_id <- paste0(string_table_1[string_table_1$search_string == x, "string_id"],
                string_table_2[string_table_2$search_string == x, "string_id"],
                combine_to[combine_to$search_string == x, "string_id"],
                collapse = ", ")
            return(collapsed_id)
        })

    collapsed_strings_table$string_id <- unlist(strings_id)
    return(collapsed_strings_table)
}

### Rewrite match string function to include both snp and gene k-mers


matched_string <- function(datum, searches){
    s_strings <- unique(searches$search_string)
    c_count <- lapply(s_strings, function(string) {
            match <- unlist(gregexpr(string, datum))
            return(match[which(match > 0)])
        })
    s_m <- lapply(c_count, length)
    matched_string <- s_strings[which(s_m > 0)]
    match_count <- unlist(s_m[which(s_m > 0)])
    return(
        data.frame(matched_string = matched_string,
            match_count = match_count,
            string_id = searches[match(matched_string, searches$search_string), "string_id"],
            stringsAsFactors = FALSE)
    )
}

combine_match_gene <- function(match_table){
    if (ncol(match_table) == 0){
        return(NULL)
    }
    combined_match_count <- aggregate(match_table$match_count,
        by=list(matched_string=match_table$matched_string), FUN=sum)
    
    strings_gene <- bplapply(match_table$string_id, function(x){
        string_gene <- strsplit(x, split = ", ")[[1]]
        string_gene <- unlist(lapply(strsplit(string_gene, split = "_"), `[`, 1))
        return(string_gene)
    })

    names(strings_gene) <- match_table$matched_string
    genes <- unique(unlist(strings_gene))

    gene_presence_indicators <- data.frame(gene = c(),
        matched_string_id = c(),
        matched_string_count = c(), stringsAsFactors = FALSE)

    for (gene in genes){
        string_has_gene <- lapply(strings_gene, function(x){
            return(gene %in% x)
        })
        gene_presence_indicators <- rbind(gene_presence_indicators, data.frame(
            gene = gene,
            matched_string_id = paste(
                names(string_has_gene[which(string_has_gene == TRUE)]), collapse = ", "),
            matched_string_count = paste(
                combined_match_count[
                    match(names(string_has_gene[which(string_has_gene == TRUE)]),
                    combined_match_count$matched_string),
                "x"], collapse = "_")
        ))
    }
    
    
    return(gene_presence_indicators)
}



library(minSNPs)
library(BiocParallel)

register(MulticoreParam(workers = 64, progressbar = TRUE))
genes <- read.csv("../gene_search_table.csv", stringsAsFactors=FALSE)

test_files <- list.files(pattern = "*.fastq")
test_data <- lapply(test_files, get_data)

final_res <- list()
i <- 1
for (dat in test_data){
    temp_res <- bplapply(dat, function(x){
        return(matched_string(paste0(x), searches = final_search))
        })
    temp_table <- do.call(rbind, temp_res)
    final_res[[i]] <- combine_match_gene(temp_table)
    i <- i + 1
}

final_search
