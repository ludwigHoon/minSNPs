setwd("/home/ludwig/MIN")
devtools::load_all("./MinSNPs")
hets_vivax <- read_fasta("n_vivax.fasta")
homs_vivax <- read_fasta("Pv_Nature_communications_2018_VQSLOD3.min0.fasta")

k2 <- c("PY0019-C", "PY0024-C", "PY0026-C", "PY0027-C", "PY0034-C",
    "PY0044-C", "PY0053-C", "PY0056-C", "PY0057-C", "PY0058-C",
    "PY0060-C", "PY0061-C", "PY0067-C", "PY0068-C", "PY0072-C",
    "PY0073-C", "PY0074-C", "PY0075-C", "PY0101-C", "PY0105-C",
    "PY0117-C", "PY0035-C", "PY0042-C", "PY0076-C", "PY0119-C", "PY0085-C")

# 10 sets of 5-SNPs with highest Simpson's score obtained from
# SNP matrix where IUPAC and failed genotyping are **not resolved**.
w_26_k2_nr_1 <- c(419073, 347252, 274608, 120632, 387620) # 0.75145910029631
w_26_k2_nr_2 <- c(488422, 347252, 274608, 120632, 387620) # 0.749902726646913
w_26_k2_nr_3 <- c(347252, 274608, 120632, 146476, 51994) # 0.571937385890874
w_26_k2_nr_4 <- c(523767, 274608, 17742, 51136, 387620) # 0.563616772919099
w_26_k2_nr_5 <- c(523988, 274608, 17742, 51136, 387620) # 0.563616772919099
w_26_k2_nr_6 <- c(341835, 274608, 17742, 51136, 387620) # 0.559995211158002
w_26_k2_nr_7 <- c(393708, 274608, 17742, 51136, 387620) # 0.559995211158002
w_26_k2_nr_8 <- c(396377, 274608, 17742, 51136, 387620) # 0.559995211158002
w_26_k2_nr_9 <- c(397652, 274608, 17742, 51136, 387620) # 0.559995211158002
w_26_k2_nr_10 <- c(397689, 274608, 17742, 51136, 387620) # 0.559995211158002



# 10 sets of 5-SNPs with highest Simpson's score obtained from
# SNP matrix where IUPAC and failed genotyping are **arbitrarily resolve**.
w_26_k2_1 <- c(3480, 182350, 466351, 168328, 286908) # 0.96947113226183
w_26_k2_2 <- c(7645, 108665, 243020, 168328, 438299) # 0.969321480949388
w_26_k2_3 <- c(27920, 77522, 312065, 168328, 259388) # 0.968603154649666
w_26_k2_4 <- c(28260, 108665, 60806, 168328, 340981) # 0.9692915506869
w_26_k2_5 <- c(36918, 178410, 370716, 168328, 154479) # 0.969171829636946
w_26_k2_6 <- c(38300, 103631, 56655, 123688, 168328) # 0.967136571787735
w_26_k2_7 <- c(39628, 184223, 257071, 168328, 140634) # 0.969501062524318
w_26_k2_8 <- c(43889, 77554, 108084, 499911, 168328) # 0.967914758612433
w_26_k2_9 <- c(47782, 267072, 493154, 161976, 168328) # 0.967615455987549
w_26_k2_10 <- c(50015, 274882, 264454, 161977, 460394) # 0.967376013887642
resolved_full <- ls(pattern = "w_26_k2_[0-9]")

# 5 sets of 5-SNPs with highest Simpson's score obtained from
# SNP matrix where IUPAC and failed genotyping are **arbitrarily resolve**,
# but only 1 randomly selected K2 isolate is in the matrix: PY0073-C.
w_1_k2_1 <- c(6476, 108505, 70773, 496877, 40623) # 0.971571108910165
w_1_k2_2 <- c(19321, 72096, 333746, 57339, 14789) # 0.971607791350281
w_1_k2_3 <- c(22213, 332870, 490690, 114468, 128937) # 0.971607791350281
w_1_k2_4 <- c(29659, 69165, 209376, 107826, 493920) # 0.97179120355086
w_1_k2_5 <- c(29798, 57336, 132649, 199609, 332975) # 0.972011298191556
resolved_1k2 <- ls(pattern = "w_1_k2_[0-9]")

# 5 sets of 5-SNPs with highest Simpson's score obtained from
# SNP matrix where IUPAC and failed genotyping are **arbitrarily resolve**,
# but only 5 randomly selected K2 isolate is in the matrix: PY0034-C, PY0061-C,
# PY0073-C, PY0042-C, PY0076-C. **Same seed as previous used**
w_5_k2_1 <- c(268, 8243, 257288, 81709, 171138) # 0.971846966634755
w_5_k2_2 <- c(35938, 129916, 415720, 335796, 291952) # 0.971917881076481
w_5_k2_3 <- c(37471, 56428, 70773, 332870, 124899) # 0.971811509413892
w_5_k2_4 <- c(41422, 56428, 75728, 9946, 37301) # 0.971705137751303
w_5_k2_5 <- c(46256, 53746, 289428, 111562, 28397) # 0.971705137751303
resolved_5k2 <- ls(pattern = "w_5_k2_[0-9]")

dir.create("./homs", showWarnings = TRUE)
for (dataset in c(resolved_full, resolved_1k2, resolved_5k2)) {
    for (data in dataset) {
        result <- find_optimised_snps(homs_vivax, metric = "percent",
            goi = k2, max_depth = 0, included_positions = get(data),
            iterate_included = TRUE)
        output_result(result, view = "csv",
            file_name = paste("./homs/", data, ".tsv", sep = ""))
    }
}
dir.create("./hets", showWarnings = TRUE)
for (dataset in c(resolved_full, resolved_1k2, resolved_5k2)) {
    for (data in dataset) {
        result <- find_optimised_snps(hets_vivax, metric = "percent",
            goi = k2, max_depth = 0, included_positions = get(data),
            iterate_included = TRUE)
        output_result(result, view = "csv",
            file_name = paste("./hets/", data, ".tsv", sep = ""))
    }
}