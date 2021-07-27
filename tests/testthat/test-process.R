
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
  expect_equal(length(names(error_file_2)), 56)
})


test_that("usual length", {
  setup()
  # Normal File
  expect_equal(usual_length(chlamydia), 19570)

  # With shorter samples, 1st and 2nd are shorter
  expect_equal(usual_length(error_file_1), 19570)
  expect_equal(usual_length(error_file_1[1]), 19569)
  expect_equal(usual_length(error_file_1[2]), 19570)
  expect_equal(usual_length(error_file_1[1:2]), 19570)
  expect_equal(usual_length(error_file_1[3]), 19570)
  expect_equal(usual_length(error_file_1[1:3]), 19570)
  expect_equal(usual_length(error_file_1[4]), 19570)
  expect_equal(usual_length(error_file_1[1:4]), 19570)

  # With both short and longer samples
  expect_equal(usual_length(error_file_2), 19570)
  expect_equal(usual_length(error_file_2[1]), 19571)
  expect_equal(usual_length(error_file_2[2]), 19570)
  expect_equal(usual_length(error_file_2[1:2]), 19571)
  expect_equal(usual_length(error_file_2[3]), 19567)
  expect_equal(usual_length(error_file_2[1:3]), 19571)
  expect_equal(usual_length(error_file_2[4]), 19570)
  expect_equal(usual_length(error_file_2[1:4]), 19570)

})
