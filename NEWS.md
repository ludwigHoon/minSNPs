# minSNPs v0.0.1
Rebuilding and cleaning old version of minSNPs.
## Changes
- `read_fasta` will be used to read fasta file, the sequence name can contain spaces.
- `write_fasta` will be used to write fasta file, the sequence name can contain spaces.
- `process_allele` can be used in place of the different functions to preprocess the FASTA file. It will perform all equivalent functions in one go; the others functions are no longer exported.
- `calculate_simpson` is used to calculate the Simpson's index from list of pattern.
- `calculate_percent` is used to calculate the dissimilarity index from list of pattern.
- `check_percent` is used to check that the necessary parameter for `calculate_percent` is passed.
- `get_metric_fun` is used internally to translate the metric to the functions, additional metric can be added by defining the functions and overriding the global `MinSNPs_metrics` variable.
- `calculate_simpson`, `calculate_percent`, and `check_percent`, `get_metric_fun` are only included as sample, there is no reason to call these functions directly.
- `find_optimised_snps` is used to find optimised SNP set(s). It can also be used to recalculate index by giving it included positions and setting max_depth to 0.