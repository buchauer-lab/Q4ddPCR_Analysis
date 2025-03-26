library(testthat)
library(dplyr)
source("../R/TotalHIVDNA.R")

# ======================== test compute_total_HIV =============================
# Sample test data frame
test_df <- tibble(data.frame(
  "intact provirus/Mio cells Pol+Gag+Psi+Env+, corrected for shearing" = c(2),
  "intact provirus/Mio cells Gag+Psi+Env+" = c(3),
  "intact provirus/Mio cells Psi+Env+" = c(5),
  "intact provirus/Mio cells Gag+Pol+Env+" = c(3),
  "Mean Target/Mio cells" = c(36),
  check.names = FALSE
))

# Test function

test_that("compute_total_HIV calculates correctly", {
  result <- compute_total_HIV(test_df)
  expect_true("total HIV DNA/Mio cells" %in% names(result))
  expected_total <- (36 / 4) - 5 + 6 - 2 # sing - doub + trip - quad
  expect_equal(result[["total HIV DNA/Mio cells"]], expected_total)
})

# Test handling missing values
test_that("compute_total_HIV handles missing values correctly", {
  df_missing <- test_df
  df_missing[1, "Mean Target/Mio cells"] <- NA
  result <- compute_total_HIV(df_missing)
  expect_true(is.na(result[["total HIV DNA/Mio cells"]]))
})

# Test when column selection results in multiple values
test_that("compute_total_HIV throws error for multiple quad values", {
  df_multiple <- test_df
  df_multiple[2, "intact provirus/Mio cells Pol+Gag+Psi+Env+shearing"] <- 7
  expect_error(
    compute_total_HIV(df_multiple),
    "More than 1 column for Total HIV calculation selected."
  )
})

# Test edge case with zero values
test_that("compute_total_HIV handles zero values correctly", {
  df_zeros <- test_df
  df_zeros[1, ] <- 0
  result <- compute_total_HIV(df_zeros)
  expect_equal(result[["total HIV DNA/Mio cells"]], 0)
})


# ======================== test compute_total_HIV_envPsi =======================
test_that("compute_total_HIV_envPsi correctly computes total HIV DNA", {
  df <- tibble(data.frame(
    Target = c("Env", "Psi", "Gag"),
    `Mean Target/Mio cells` = c(10, 15, 5),
    `intact provirus/Mio cells Psi.Env, corrected for shearing` = c(7, NA, 7),
    `intact provirus/Mio cells Env.Psi, corrected for shearing` = c(NA, NA, NA), # Fallback column empty
    check.names = FALSE
  ))
  
  df_result <- compute_total_HIV_envPsi(df)
  
  # Expected total HIV DNA = env + psi - envpsi = 10 + 15 - 7 = 18
  expect_equal(unique(df_result[["total HIV DNA/Mio cells (Env.Psi)"]]), 18)
})

test_that("compute_total_HIV_envPsi handles missing EnvPsi values correctly", {
  df <- data.frame(
    Target = c("Env", "Psi", "Gag"),
    `Mean Target/Mio cells` = c(20, 25, 15),
    #`intact provirus/Mio cells Psi.Env, corrected for shearing` = c(NA, NA, NA), # Empty
    `intact provirus/Mio cells Env.Psi, corrected for shearing` = c(12, NA, 12), # Uses fallback
    check.names = FALSE
  )
  
  df_result <- compute_total_HIV_envPsi(df)
  
  # Expected total HIV DNA = env + psi - envpsi = 20 + 25 - 12 = 33
  expect_equal(unique(df_result[["total HIV DNA/Mio cells (Env.Psi)"]]), 33)
})

test_that("compute_total_HIV_envPsi correctly handles cases where both EnvPsi columns are missing", {
  df <- data.frame(
    Target = c("Env", "Psi"),
    `Mean Target/Mio cells` = c(5, 10),
    #`intact provirus/Mio cells Psi.Env, corrected for shearing` = c(NA, NA),
    #`intact provirus/Mio cells Env.Psi, corrected for shearing` = c(NA, NA),
    check.names = FALSE
  )
  
  expect_warning(df_result <- compute_total_HIV_envPsi(df),
                 "No Env+Psi+ detected. Will use 0.",
                 fixed=TRUE)
  
  # Expected total HIV DNA = env + psi - envpsi = 5 + 10 - 0 = 15 (envpsi defaults to 0)
  expect_equal(unique(df_result[["total HIV DNA/Mio cells (Env.Psi)"]]), 15)
})

test_that("compute_total_HIV_envPsi throws an error for missing columns", {
  df <- data.frame(
    Target = c("Env", "Psi"),
    `intact provirus/Mio cells Psi.Env, corrected for shearing` = c(5, NA),
    check.names = FALSE
  )
  
  expect_error(compute_total_HIV_envPsi(df), "No Mean Target/Mio cells detected")
})

test_that("compute_total_HIV_envPsi handles an empty dataframe gracefully", {
  df <- data.frame(Target = character(), `Mean Target/Mio cells` = numeric(),
                   `intact provirus/Mio cells Psi.Env, corrected for shearing` = numeric(),
                   check.names = FALSE)
  
  df_result <- compute_total_HIV_envPsi(df)
  
  expect_true(nrow(df_result) == 0)
  expect_true("total HIV DNA/Mio cells (Env.Psi)" %in% colnames(df_result))
})

