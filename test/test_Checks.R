#test functions in readers script
source("../R/Checks.R")
library(testthat)


# =============== test is_named_vector =========================================
test_that("is_named_vector correctly identifies named vectors", {
  expect_true(is_named_vector(c(a = 1, b = 2, c = 3)))  # Named numeric vector
  expect_true(is_named_vector(setNames(c("x", "y"), c("name1", "name2"))))  # Named character vector
})

test_that("is_named_vector rejects unnamed vectors", {
  expect_false(is_named_vector(c(1, 2, 3)))  # Unnamed numeric vector
  expect_false(is_named_vector(c("x", "y")))  # Unnamed character vector
})

test_that("is_named_vector rejects non-atomic objects", {
  expect_false(is_named_vector(list(a = 1, b = 2)))  # Named list
  expect_false(is_named_vector(data.frame(a = 1, b = 2)))  # Data frame
})

test_that("is_named_vector rejects NULL and empty values", {
  expect_false(is_named_vector(NULL))  # NULL value
  expect_false(is_named_vector(c()))  # Empty vector
})

test_that("is_named_vector rejects vectors with empty names", {
  bad_vec <- c(1, 2, 3)
  names(bad_vec) <- c("a", "", "c")  # One missing name
  expect_false(is_named_vector(bad_vec))
})
