#test functions in readers script
source("./R/CreateHouseholdTable.R")
library(testthat)

# =============== test define_groups ===========================================
test_that("define_groups correctly groups wells", {
  dtQC <- data.frame(
    Well = c("A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample1", "Sample1", "Sample2", "Sample2", "Sample2", "Sample2"),
    Target = c("RPP30", "RPP30Shear", "RPP30", "RPP30Shear", "RPP30", "RPP30Shear", "RPP30", "RPP30Shear"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  dilution_factor <- setNames(c(2, 3), c("Sample1", "Sample2"))
  
  result <- define_groups(dtQC, dilution_factor)
  
  expect_true("shared_Wells" %in% colnames(result))
  expect_equal(nrow(result), 2)
  expect_equal(result$shared_Wells[1], "A1, A2")
  expect_equal(result$shared_Wells[2], "B1, B2")
  expect_equal(unname(result$`Dilution Factor`[1]), 2)
  expect_equal(unname(result$`Dilution Factor`[2]), 3)
})

test_that("define_groups throws an error when RPP30 and RPP30Shear do not match", {
  dtQC <- data.frame(
    Well = c("A1", "A2", "B1"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample2"),
    Target = c("RPP30", "RPP30Shear", "RPP30"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  expect_error(define_groups(dtQC, NULL), "Wells with targets RPP30 and RPP30Shear do not match")
})

test_that("define_groups warns when dilution factor is missing", {
  dtQC <- data.frame(
    Well = c("A1", "A1", "A2", "A2"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample1", "Sample1"),
    Target = c("RPP30", "RPP30Shear", "RPP30", "RPP30Shear"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  expect_warning(result <- define_groups(dtQC, NULL), "No dilution factor has been specified")
  expect_equal(result$`Dilution Factor`, rep(1, nrow(result)))
})


# =============== test add_dilution_factor =====================================
test_that("add_dilution_factor correctly assigns dilution factors", {
  df <- data.frame(
    `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  dilution_factor <- setNames(c(2, 3, 4), c("Sample1", "Sample2", "Sample3"))
  
  result <- add_dilution_factor(df, dilution_factor)
  
  expect_equal(result$`Dilution Factor`, c(2, 3, 4))
})

test_that("add_dilution_factor assigns default dilution factor of 1 when missing", {
  df <- data.frame(
    `Sample description 1` = c("Sample1", "Sample2", "Sample3"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  dilution_factor <- setNames(c(2), c("Sample1"))
  
  expect_warning(result <- add_dilution_factor(df, dilution_factor), "Assuming 1 as dilution factor if not specified otherwise.")
  expect_equal(result$`Dilution Factor`, c(2, 1, 1))
})

test_that("add_dilution_factor warns when dilution factor names additionally do not match sample descriptions", {
  df <- data.frame(
    `Sample description 1` = c("Sample1", "Sample3"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  dilution_factor <- setNames(c(2, 3, 4), c("Sample1", "Sample2", "SampleX"))
  
  expect_warning(add_dilution_factor(df, dilution_factor), "Not all specified dilution factor names appear in sample description.")
})

# =============== test sufficient_droplets =====================================
test_that("sufficient_droplets correctly removes wells below the threshold", {
  df <- data.frame(
    Well = c("A1", "A2", "B1", "B2"),
    `Accepted Droplets` = c(10, 5, 12, 3),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  expect_warning(result <- sufficient_droplets(df, thresh = 6), "Removed well")
  
  expect_equal(nrow(result), 2)  # Two rows should remain (A1, B1)
  expect_true(all(result$`Accepted Droplets` >= 6))  # All remaining rows should have droplets >= 6
})

test_that("sufficient_droplets issues a warning when removing wells below the threshold", {
  df <- data.frame(
    Well = c("A1", "A2", "B1", "B2"),
    `Accepted Droplets` = c(10, 5, 12, 3),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  expect_warning(
    result <- sufficient_droplets(df, thresh = 6),
    "Removed well(s): A2, B2: Less than 6 droplets.",
    fixed = TRUE
  )
})

test_that("sufficient_droplets keeps all wells if they meet the threshold", {
  df <- data.frame(
    Well = c("A1", "A2", "B1", "B2"),
    `Accepted Droplets` = c(10, 15, 12, 20),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  result <- sufficient_droplets(df, thresh = 6)
  
  expect_equal(nrow(result), 4)  # All rows should remain as they meet the threshold
  expect_true(all(result$`Accepted Droplets` >= 6))  # All rows should have droplets >= 6
})

test_that("sufficient_droplets works with custom droplet column", {
  df <- data.frame(
    Well = c("A1", "A2", "B1", "B2"),
    Custom_Droplets = c(10, 5, 12, 3),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  expect_warning(result <- sufficient_droplets(df, thresh = 6, droplet_col = "Custom_Droplets"),
                 "Removed well(s): A2, B2: Less than 6 droplets.",
                 fixed = TRUE)
  
  expect_equal(nrow(result), 2)  # Two rows should remain (A1, B1)
  expect_true(all(result$Custom_Droplets >= 6))  # All remaining rows should have droplets >= 6
})

# =============== test create_household_table ==================================
test_that("create_household_table correctly computes the table and grouped data", {
  # Sample input data
  dtQC <- data.frame(
    Well = c("A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample1", "Sample1",
                               "Sample2", "Sample2", "Sample2", "Sample2"),
    Target = c("RPP30", "RPP30Shear", "RPP30", "RPP30Shear", 
               "RPP30", "RPP30Shear", "RPP30", "RPP30Shear"),
    `Conc(copies/µL)` = c(100, 150, 200, 250, 100, 150, 200, 250),
    `Ch1+Ch2+` = c(1000, 1200, 800, 900, 1000, 1200, 800, 900),
    `Ch1+Ch2-` = c(300, 320, 250, 270, 300, 320, 250, 270),
    `Ch1-Ch2+` = c(200, 220, 180, 190, 200, 220, 180, 190),
    `Accepted Droplets` = c(10, 10, 10, 10, 10, 10, 10, 10),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  grouped_data <- data.frame(
    `Sample description 1` = c("Sample1", "Sample2"),
    `Dilution Factor` = c(2, 3),
    RPP30 = c(TRUE, TRUE),
    shared_Wells = c("A1, A2", "B1, B2"),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  # Parameters
  thresh <- 5
  mean_copies_factor <- 1.5
  mean_cells_per_reac_factor <- 0.5
  
  # Run the function
  result <- create_household_table(dtQC, grouped_data, thresh, mean_copies_factor, mean_cells_per_reac_factor)
  
  # Extract results
  tab1 <- result[[1]]
  grouped_data_updated <- result[[2]]
  
  # Tests
  expect_equal(nrow(tab1), 8)  # All 4 rows should be included as droplets > thresh
  expect_true("Mean concentration RPP30 + RPP30Shear (copies/µL)" %in% colnames(tab1))  # Check if mean concentration is calculated
  expect_true("Shearing index" %in% colnames(tab1))  # Check if shearing index is calculated
  expect_true("Mean cells per reaction" %in% colnames(tab1))  # Check if cells per reaction is calculated
  
  # Check that grouped_data is updated with the expected columns
  expect_true("Mean concentration RPP30 + RPP30Shear (copies/µL)" %in% colnames(grouped_data_updated))
  expect_true("Mean unsheared" %in% colnames(grouped_data_updated))
  expect_true("Mean cells per reaction" %in% colnames(grouped_data_updated))
})

test_that("create_household_table removes wells with insufficient droplets", {
  # Sample input data
  dtQC <- data.frame(
    Well = c("A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample1", "Sample1",
                               "Sample2", "Sample2", "Sample2", "Sample2"),
    Target = c("RPP30", "RPP30Shear", "RPP30", "RPP30Shear", 
               "RPP30", "RPP30Shear", "RPP30", "RPP30Shear"),
    `Conc(copies/µL)` = c(100, 150, 200, 250, 100, 150, 200, 250),
    `Ch1+Ch2+` = c(1000, 1200, 800, 900, 1000, 1200, 800, 900),
    `Ch1+Ch2-` = c(300, 320, 250, 270, 300, 320, 250, 270),
    `Ch1-Ch2+` = c(200, 220, 180, 190, 200, 220, 180, 190),
    `Accepted Droplets` = c(10, 10, 10, 10, 10, 10, 9, 9),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  grouped_data <- data.frame(
    `Sample description 1` = c("Sample1", "Sample2"),
    `Dilution Factor` = c(2, 3),
    RPP30 = c(TRUE, TRUE),
    shared_Wells = c("A1, A2", "B1, B2"),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  # Set threshold to 10 droplets, so B2 should be removed
  thresh <- 10
  mean_copies_factor <- 1.5
  mean_cells_per_reac_factor <- 0.5
  
  # Run the function
  expect_warning(result <- create_household_table(dtQC, grouped_data, thresh,
                                                  mean_copies_factor, mean_cells_per_reac_factor),
                 "Removed well(s): B2: Less than 10 droplets.",
                 fixed=TRUE)
  
  # Extract results
  tab1 <- result[[1]]
  
  # Test that B2 is removed due to droplets being below the threshold
  expect_equal(nrow(tab1), 6)  # Only 3 wells should remain
  expect_true(!"B2" %in% tab1$Well)  # Check B2 is removed
})


test_that("create_household_table computes correct shearing index", {
  # Sample input data
  dtQC <- data.frame(
    Well = c("A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"),
    `Sample description 1` = c("Sample1", "Sample1", "Sample1", "Sample1",
                               "Sample2", "Sample2", "Sample2", "Sample2"),
    Target = c("RPP30", "RPP30Shear", "RPP30", "RPP30Shear", 
               "RPP30", "RPP30Shear", "RPP30", "RPP30Shear"),
    `Conc(copies/µL)` = c(100, 150, 200, 250, 100, 150, 200, 250),
    `Ch1+Ch2+` = c(1000, 1200, 800, 900, 1000, 1200, 800, 900),
    `Ch1+Ch2-` = c(300, 320, 250, 270, 300, 320, 250, 270),
    `Ch1-Ch2+` = c(200, 220, 180, 190, 200, 220, 180, 190),
    `Accepted Droplets` = c(10, 10, 10, 10, 10, 10, 9, 9),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  grouped_data <- data.frame(
    `Sample description 1` = c("Sample1", "Sample2"),
    `Dilution Factor` = c(2, 3),
    RPP30 = c(TRUE, TRUE),
    shared_Wells = c("A1, A2", "B1, B2"),
    check.names = FALSE,  # Added check.names = FALSE
    stringsAsFactors = FALSE
  )
  
  # Set threshold to 5
  thresh <- 5
  mean_copies_factor <- 1.5
  mean_cells_per_reac_factor <- 0.5
  
  # Run the function
  result <- create_household_table(dtQC, grouped_data, thresh, mean_copies_factor, mean_cells_per_reac_factor)
  
  # Extract results
  tab1 <- result[[1]]
  
  # Test if Shearing index is correctly calculated
  expect_true("Shearing index" %in% colnames(tab1))  # Check if Shearing index exists
  expect_true(all(tab1$`Shearing index` >= 0))  # Shearing index should be between 0 and 1
})

















