
setup <- function() {
  chlamydia <<- read.fasta(file =
    system.file("extdata", "Chlamydia_mapped.fasta", package = "minSNPs"))
  error_file_1 <<- read.fasta(file =
    system.file("extdata", "Chlamydia_1.fasta", package = "minSNPs"))
  error_file_2 <<- read.fasta(file =
    system.file("extdata", "Chlamydia_2.fasta", package = "minSNPs"))
}

test_that("read fasta file", {
  setup()
  expect_equal(length(names(chlamydia)), 56)
  expect_equal(length(names(error_file_1)), 56)
  expect_equal(length(names(error_file_2)), 58)
})


test_that("usual length", {
  setup()
  # Normal File
  expect_equal(get_usual_length(chlamydia), 19570)

  # With shorter samples, 1st and 11th are shorter
  expect_equal(get_usual_length(error_file_1), 19570)
  expect_equal(get_usual_length(error_file_1[1]), 19569)
  expect_equal(get_usual_length(error_file_1[2]), 19570)
  expect_equal(get_usual_length(error_file_1[1:2]), 19570)
  expect_equal(get_usual_length(error_file_1[3]), 19570)
  expect_equal(get_usual_length(error_file_1[1:3]), 19570)
  expect_equal(get_usual_length(error_file_1[4]), 19570)
  expect_equal(get_usual_length(error_file_1[1:4]), 19570)
  expect_equal(get_usual_length(error_file_1[10]), 19570)
  expect_equal(get_usual_length(error_file_1[11]), 19569)
  expect_equal(get_usual_length(error_file_1[10:11]), 19570)
  expect_equal(get_usual_length(error_file_1[12]), 19570)
  expect_equal(get_usual_length(error_file_1[11:12]), 19570)
  expect_equal(get_usual_length(error_file_1[10:12]), 19570)

  # With both short and longer samples
  expect_equal(get_usual_length(error_file_2), 19570)
  expect_equal(get_usual_length(error_file_2[1]), 19571)
  expect_equal(get_usual_length(error_file_2[2]), 19570)
  expect_equal(get_usual_length(error_file_2[1:2]), 19571)
  expect_equal(get_usual_length(error_file_2[3]), 19567)
  expect_equal(get_usual_length(error_file_2[1:3]), 19571)
  expect_equal(get_usual_length(error_file_2[4]), 19570)
  expect_equal(get_usual_length(error_file_2[1:4]), 19570)
  expect_equal(get_usual_length(error_file_2[55]), 19570)
  expect_equal(get_usual_length(error_file_2[56]), 19563)
  expect_equal(get_usual_length(error_file_2[55:56]), 19570)
})

test_that("flagging allele", {
  setup()
  expect_vector(flag_allele(chlamydia), ptype = character(), size = 0)

  expect_vector(flag_allele(error_file_1), ptype = character(), size = 2)
  expect_vector(flag_allele(error_file_1[1:10]), ptype = character(), size = 1)
  expect_vector(flag_allele(error_file_1[11:20]), ptype = character(), size = 1)
  expect_equal(flag_allele(error_file_1), c("A_D213", "Ba_Aus25"))

  expect_vector(flag_allele(error_file_2), ptype = character(), size = 5)
  expect_vector(flag_allele(error_file_2[1:10]), ptype = character(), size = 2)
  expect_vector(flag_allele(error_file_2[50:56]), ptype = character(), size = 1)
  expect_equal(unique(flag_allele(error_file_2)),
    c("A_D213", "Ia_SotoGIa3", "D_SotoGD1"))
})

test_that("flagging duplicated allele name", {
  setup()
  expect_equal(remove_dup_allele(chlamydia), chlamydia)
  expect_equal(remove_dup_allele(error_file_1), error_file_1)
  expect_equal(length(remove_dup_allele(error_file_2)), 56)
  expect_equal(remove_dup_allele(error_file_2)["A_D213"], error_file_2[1])
})

test_that("flagging position", {
  setup()
  expect_vector(flag_position(chlamydia), ptype = integer(), size = 0)
  # By default read.fasta force to lower case
  expect_vector(flag_position(chlamydia, ignore_case = FALSE),
    ptype = integer(), size = 19570)
  expect_vector(flag_position(chlamydia, accepted_char = c("A", "C", "T", "G"),
    ignore_case = FALSE), ptype = integer(), size = 19570)
  expect_equal(flag_position(chlamydia[1], accepted_char = c("A", "C", "T"),
    ignore_case = TRUE), which(chlamydia[[1]] == "g"))

  expect_equal(flag_position(error_file_1), c(22, 24))
  expect_equal(flag_position(error_file_1, dash_ignore = FALSE), c(22))

  expect_equal(flag_position(error_file_2[1:56]), c(1, 2, 3, 4, 120, 121))
  expect_equal(flag_position(error_file_2[1:56], dash_ignore = FALSE),
    c(1, 2, 120, 121))
})

test_that("Overall", {
  setup()
  expect_equal(length(process_allele(chlamydia)[["seqc"]]), 56)
  expect_vector(process_allele(chlamydia)[["ignored_allele"]],
    ptype = character(), size = 0)
  expect_vector(process_allele(chlamydia)[["ignored_position"]],
    ptype = integer(), size = 0)

  expect_equal(length(process_allele(error_file_1)[["seqc"]]), 54)
  expect_vector(process_allele(error_file_1)[["ignored_allele"]],
    ptype = character(), size = 2)
  expect_vector(process_allele(error_file_1)[["ignored_position"]],
    ptype = integer(), size = 2)

  expect_equal(length(process_allele(error_file_2)[["seqc"]]), 53)
  expect_vector(process_allele(error_file_2)[["ignored_allele"]],
    ptype = character(), size = 3)
  expect_true(all(process_allele(error_file_2)[["ignored_position"]] %in%
    c(1, 2, 3, 4, 120, 121)))
})