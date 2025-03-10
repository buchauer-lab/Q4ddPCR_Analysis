#test functions in readers script
source("../R/Readers.R")
library(testthat)

# =============== test read_xlsx ==============================================
test_that("test read_xlsx", {
  expect_error(read_xlsx("README.md"))
  expect_error(read_xlsx("NotExistingFile.xlsx"), 
                          regexp = "Specified xlsx file NotExistingFile.xlsx does not exist. Please check if the path and name are correct and that there are no spelling mistakes.")
  expect_type(read_xlsx("../../data/TestData/DataSheet19.06.xlsx"), "list")
  expect_visible(read_xlsx("../../data/TestData/DataSheet19.06.xlsx"))
  expect_s3_class(read_xlsx("../../data/TestData/DataSheet19.06.xlsx"),"data.frame")
})

# =============== test read_csv ===============================================
test_that("test read_csv", {
  expect_error(read_csv("README.md", 3))
  expect_error(read_csv("NotExistingFile.csv", 3),
                         regexp = paste0("Specified csv file NotExistingFile.csv does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  expect_type(read_csv("../../data/TestData/ClusterData19.06.csv", 4), "list")
  expect_visible(read_csv("../../data/TestData/ClusterData19.06.csv", 4))
  expect_error(read_csv("../../data/TestData/ClusterData19.06.csv", 2))
  expect_s3_class(read_csv("../../data/TestData/ClusterData19.06.csv", 4),"data.frame")
})

# =============== test rm_channel ==============================================

test_that("rm_channels removes specified wells", {
  df <- data.frame(Well = c("A1", "B2", "C3", "D4"), Value = c(10, 20, 30, 40))
  channels_to_remove <- c("A1", "C3")
  result <- rm_channels(df, channels_to_remove)
  
  expected_df <- data.frame(Well = c("B2", "D4"), Value = c(20, 40))
  expect_true(all(result == expected_df))
})

test_that("rm_channels returns the same dataframe when no channels match", {
  df <- data.frame(Well = c("A1", "B2", "C3"), Value = c(10, 20, 30))
  channels_to_remove <- c("X9")
  expect_error(rm_channels(df, channels_to_remove),
               "At least one of the Wells that should be removed does not occur")
})

test_that("rm_channels stops when Well column is missing", {
  df <- data.frame(ID = c("A1", "B2", "C3"), Value = c(10, 20, 30))
  channels_to_remove <- c("A1")
  expect_error(rm_channels(df, channels_to_remove),
               "No Wells specified in data, thus cannot be removed")
})

test_that("rm_channels removes all channels when all match", {
  df <- data.frame(Well = c("A1", "B2"), Value = c(10, 20))
  channels_to_remove <- c("A1", "B2")
  result <- rm_channels(df, channels_to_remove)
  
  expected_df <- data.frame(Well = character(0), Value = numeric(0))
  expect_equal(result, expected_df)
})

# =============== test create_dtQC =============================================

test_that("create_dtQC extracts correct columns", {
  df <- data.frame(`Well` = c("A1", "B2"), 
                   `Sample description 1` = c("Sample1", "Sample2"),
                   `DyeName(s)` = c("FAM", "VIC"),
                   `Target` = c("Gene1", "Gene2"),
                   `Conc(copies/µL)` = c(100, 200),
                   `Accepted Droplets` = c(3000, 6000),
                   `Positives` = c(150, 300),
                   `Negatives` = c(2850, 5700),
                   check.names = FALSE)
  
  result <- create_dtQC(df)
  expect_true(all(c("Well", "Sample description 1", "DyeName(s)", "Target", 
                    "Conc(copies/µL)", "Accepted Droplets", "Positives", 
                    "Negatives", "Threshold", "Total positives") %in% names(result)))
})

test_that("create_dtQC correctly calculates threshold", {
  df <- data.frame(Well = c("A1"), `Sample description 1` = c("Sample1"), `DyeName(s)` = c("FAM"), 
                   Target = c("Gene1"), `Conc(copies/µL)` = c(100), `Accepted Droplets` = c(3000),
                   Positives = c(150), Negatives = c(2850),
                   check.names = FALSE)
  result <- create_dtQC(df)
  expect_equal(result$Threshold, 1000)
})

test_that("create_dtQC correctly calculates total positives", {
  df <- data.frame(Well = c("A1"), `Sample description 1` = c("Sample1"), `DyeName(s)` = c("FAM"), 
                   Target = c("Gene1"), `Conc(copies/µL)` = c(100), `Accepted Droplets` = c(3000),
                   Positives = c(150), Negatives = c(2850), `Ch1+Ch2+` = c(50), `Ch1-Ch2+` = c(100),
                   `Ch1-Ch2` = c(20),
                   check.names = FALSE)
  result <- create_dtQC(df)
  expect_equal(result$`Total positives`, 150)
})

test_that("create_dtQC stops when required columns are missing", {
  df <- data.frame(Well = c("A1"), `Accepted Droplets` = c(3000))
  expect_error(create_dtQC(df), "Some required columns are missing from the data frame.")
})

# =============== test rm_zero_channel =========================================
test_that("rm_zero_channel removes wells with zero concentration", {
  dtQC <- data.frame(Well = c("A1", "B2", "C3"), 
                     `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
                     `Conc(copies/µL)` = c(100, 0, 200),
                     check.names = FALSE)
  in_csv <- data.frame(Well = c("A1", "B2", "C3"),
                       `Target 1`= c(0, 1, 2),
                       check.names = FALSE)
  
  expect_warning(result <- rm_zero_channel(dtQC, in_csv), "These wells will be removed")
  expect_false("B2" %in% result[[1]]$Well)
  expect_false("B2" %in% result[[2]]$Well)
})

test_that("rm_zero_channel retains H2O wells", {
  dtQC <- data.frame(Well = c("A1", "B2", "C3"), 
                     `Sample description 1` = c("Sample1", "H2O", "Sample3"),
                     `Conc(copies/µL)` = c(100, 0, 200),
                     check.names = FALSE)
  in_csv <- data.frame(Well = c("A1", "B2", "C3"),
                       `Target 1`= c(0, 1, 2),
                       check.names = FALSE)
  
  result <- rm_zero_channel(dtQC, in_csv)
  expect_true("B2" %in% result[[1]]$Well)
  expect_true("B2" %in% result[[2]]$Well)
})

test_that("rm_zero_channel stops when required columns are missing", {
  dtQC <- data.frame(Well = c("A1"), `Accepted Droplets` = c(3000), check.names = FALSE)
  in_csv <- data.frame(Well = c("A1"))
  expect_error(rm_zero_channel(dtQC, in_csv), "No concentration column ")
})

# =============== test read_files ==============================================

test_that("read_files reads and processes data correctly", {
  csv_data <- data.frame(Well = c("A1", "B2", "C3", ""), Value = c(10, 20, 30, ""),
                         empty = c(NA, NA, NA, NA))
  xlsx_data <- data.frame(Well = c("A1", "B2", "C3"),
                          `Conc(copies/µL)` = c(100, 200, 300),
                          `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
                          `DyeName(s)` = c("x", "y", "z"),
                          `Target` = c("I", "am", "Groot"),
                          `Accepted Droplets` = c(40, 40, 20),
                          `Positives` = c(39, 38, 19),
                          `Negatives` = c(1, 2, 1),
                          check.names = FALSE)
  
  write.csv(csv_data, "test.csv", row.names = FALSE)
  openxlsx::write.xlsx(xlsx_data, "test.xlsx")
  
  result <- read_files("test.xlsx", "test.csv", 0)
  
  expect_true(is.list(result))
  expect_true(all(c("Well", "Value") %in% names(result[[1]])))
  expect_true(all(c("Well", "Conc(copies/µL)", "Sample description 1") %in% names(result[[2]])))
})

test_that("read_files removes specified channels", {
  csv_data <- data.frame(Well = c("A1", "B2", "C3", ""), Value = c(10, 20, 30, ""),
                         empty = c(NA, NA, NA, NA))
  xlsx_data <- data.frame(Well = c("A1", "B2", "C3"),
                          `Conc(copies/µL)` = c(100, 200, 300),
                          `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
                          `DyeName(s)` = c("x", "y", "z"),
                          `Target` = c("I", "am", "Groot"),
                          `Accepted Droplets` = c(40, 40, 20),
                          `Positives` = c(39, 38, 19),
                          `Negatives` = c(1, 2, 1),
                          check.names = FALSE)
  
  write.csv(csv_data, "test.csv", row.names = FALSE)
  openxlsx::write.xlsx(xlsx_data, "test.xlsx")
  
  result <- read_files("test.xlsx", "test.csv", 0, remove_channel = c("B2"))
  
  expect_false("B2" %in% result[[1]]$Well)
  expect_false("B2" %in% result[[2]]$Well)
})

test_that("read_files removes zero concentration channels when specified", {
  csv_data <- data.frame(Well = c("A1", "B2", "C3", ""), Value = c(10, 20, 30, ""),
                         empty = c(NA, NA, NA, NA))
  xlsx_data <- data.frame(Well = c("A1", "B2", "C3"),
                          `Conc(copies/µL)` = c(100, 0, 300),
                          `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
                          `DyeName(s)` = c("x", "y", "z"),
                          `Target` = c("I", "am", "Groot"),
                          `Accepted Droplets` = c(40, 40, 20),
                          `Positives` = c(39, 38, 19),
                          `Negatives` = c(1, 2, 1),
                          check.names = FALSE)
  
  write.csv(csv_data, "test.csv", row.names = FALSE)
  openxlsx::write.xlsx(xlsx_data, "test.xlsx")

  
  expect_warning(result <- read_files("test.xlsx", "test.csv", 0, rm_zero_channels = TRUE),
                 "These wells will be removed,")
  
  expect_false("B2" %in% result[[1]]$Well)
  expect_false("B2" %in% result[[2]]$Well)
})


# =============== test read_multiple_files =====================================
test_that("read_multiple_files processes multiple files correctly", {
  csv_files <- c("test1.csv", "test2.csv")
  xlsx_files <- c("test1.xlsx", "test2.xlsx")
  
  for (i in 1:2) {
    csv_data <- data.frame(Well = c("A1", "B2", "C3", ""), Value = c(10, 20, 30, ""),
                           empty = c(NA, NA, NA, NA))
    xlsx_data <- data.frame(Well = c("A1", "B2", "C3"),
                            `Conc(copies/µL)` = c(100, 0, 300),
                            `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
                            `DyeName(s)` = c("x", "y", "z"),
                            `Target` = c("I", "am", "Groot"),
                            `Accepted Droplets` = c(40, 40, 20),
                            `Positives` = c(39, 38, 19),
                            `Negatives` = c(1, 2, 1),
                            check.names = FALSE)
    
    write.csv(csv_data, csv_files[i], row.names = FALSE)
    openxlsx::write.xlsx(xlsx_data, xlsx_files[i])
  }
  
  result <- read_multiple_files(xlsx_files, csv_files, csv_skip = 0)
  
  expect_true(is.list(result))
  expect_equal(nrow(result[[1]]), 6)
  expect_equal(nrow(result[[2]]), 6)
})












