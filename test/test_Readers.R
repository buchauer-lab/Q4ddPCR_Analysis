#test functions in readers script
source("R/Readers.R")
library(testthat)

# =============== test read_xlsx ==============================================
test_that("test read_xlsx", {
  expect_error(read_xlsx("README.md"))
  expect_error(read_xlsx("NotExistingFile.xlsx"), 
                          regexp = "Specified xlsx file NotExistingFile.xlsx does not exist. Please check if the path and name are correct and that there are no spelling mistakes.")
  expect_type(read_xlsx("../data/TestData/DataSheet19.06.xlsx"), "list")
  expect_visible(read_xlsx("../data/TestData/DataSheet19.06.xlsx"))
  expect_s3_class(read_xlsx("../data/TestData/DataSheet19.06.xlsx"),"data.frame")
})

# =============== test read_csv ===============================================
test_that("test read_csv", {
  expect_error(read_csv("README.md", 3, 285))
  expect_error(read_csv("NotExistingFile.csv", 3, 285),
                         regexp = paste0("Specified csv file NotExistingFile.csv does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  expect_type(read_csv("../data/TestData/ClusterData19.06.csv", 4, 285), "list")
  expect_visible(read_csv("../data/TestData/ClusterData19.06.csv", 4, 285))
  expect_error(read_csv("../data/TestData/ClusterData19.06.csv", 2, 285))
  expect_error(read_csv("../data/TestData/ClusterData19.06.csv", 4, 288), 
               regexp = "Too many rows from csv are read. Potentially 3 too many.")
  expect_s3_class(read_csv("../data/TestData/ClusterData19.06.csv", 4, 285),"data.frame")
})
