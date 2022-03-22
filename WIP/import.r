process_result_file <- function(result_file) {
    con <- file(result_file, "r")

    new_result <- F
    result <- data.frame(SNPs = c(), result = c())

    while(TRUE){
        data <- readLines(con, n = 1)
        if (length(data) == 0){
            break
        }

        if (! new_result) {
            if (startsWith(data, "Result")) {
                new_result <- T
                data <- readLines(con, n = 1)
            }            
        }
        if (new_result) {
            while (TRUE) {
                data <- readLines(con, n = 1)
                if (data != ""){
                    prev <- data
                } else {
                    break
                }
            }
            splitted <- strsplit(prev, split = "\t")[[1]]
            res <- as.numeric(splitted[2])
            snps <- splitted[1]
            result <- rbind(result, data.frame(SNPs = snps, result = res))
            new_result <- F
        }
    }
    close(con)
    return(result)
}