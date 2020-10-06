#https://kbroman.org/pkg_primer/pages/tests.html
#https://r-pkgs.org/tests.html
# TO DO: context("Checking Input files")

context("String length")


test_that("str_length is number of characters", {
  expect_equal(stringr::str_length("a"), 1)
  expect_equal(stringr::str_length("ab"), 2)
  expect_equal(stringr::str_length("abc"), 3)
})


test_that("str_length of factor is length of level", {
  expect_equal(stringr::str_length(factor("a")), 1)
  expect_equal(stringr::str_length(factor("ab")), 2)
  expect_equal(stringr::str_length(factor("abc")), 3)
})


test_that("str_length of missing is missing", {
  expect_equal(stringr::str_length(NA), NA_integer_)
  expect_equal(stringr::str_length(c(NA, 1)), c(NA, 1))
  expect_equal(stringr::str_length("NA"), 2)
})

