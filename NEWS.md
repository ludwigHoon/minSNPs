# minSNPs v0.1.0
- Added method to generate search strings (`generate_snp_search_string`) from identified SNPs set.
- Added method to search through long read sequence (`search_from_fastq_reads`).

# minSNPs v0.0.5
Fixed bug in `find_optimised_snps` where when multiple SNPs sets are looked, not all SNPs set will iterate to the next depth even if max_depth has not been reached and resolving power is not 1.

# minSNPs v0.0.4
Fixed bug in `find_optimised_snps` when there are included SNPs.

# minSNPs v0.0.3
Minor changes & bug removal

## Changes
- `find_optimised_snps` can use `BiocParallel::SnowParam()` for windows.

# minSNPs v0.0.2
Minor changes & testing

## Changes
- `find_optimised_snps` can receive an additional argument `output_progress = TRUE`, such that progress can be shown in terminal when looking for multiple result.
- `write_fasta` bug where only the first letter is written to the file is fixed.
- `output_result` bug where residuals are not returned.
- `process_allele` will print names of duplicated isolate, and keep only the first instance of the found isolate.
- `iterate_merge` can be used to merge alignment matrices.
- `find_optimised_snps` reduced memory usage.
- `resolve_IUPAC_missing` can be used to substitute the ambiguity codes found in the sequence. 
- `output_result` updated such that it can be presented correctly in excel.
- `cal_fn` and `cal_fp` can be used to calculate and identify false negatives and false postives respectively. To be called with `cal_met_snp`

# minSNPs v0.0.1
First release to CRAN.

## Changes
- `read_fasta` will be used to read fasta file, the sequence name can contain spaces.
- `write_fasta` will be used to write fasta file, the sequence name can contain spaces.
- `process_allele` can be used in place of the different functions to preprocess the FASTA file. It will perform all equivalent functions in one go; the others functions are no longer exported.
- `calculate_simpson` is used to calculate the Simpson's index from list of pattern.
- `calculate_percent` is used to calculate the dissimilarity index from list of pattern.
- `check_percent` is used to check that the necessary parameter for `calculate_percent` is passed.
- `get_metric_fun` is used to translate the metric to the functions, additional metric can be added by defining the functions and overriding the global `MinSNPs_metrics` variable.
- `calculate_simpson`, `calculate_percent`, and `check_percent`, `get_metric_fun` are only as default, there is no reason to call these functions directly.
- `find_optimised_snps` is used to find optimised SNP set(s). It can also be used to recalculate index by giving it included positions and setting max_depth to 0.
- `output_result` is used to output results from `find_optimised_snps`.