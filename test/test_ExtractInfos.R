library(testthat)
library(dplyr)
library(mockery)
library(tidyr)
source("../R/ExtractInfos.R")

# =============== test get_multipos ===========================================
test_that("get_multipos works correctly", {
  result <- get_multipos(c("A", "B", "C"))

  expect_equal(length(result), 4)
  expect_equal(unlist(result), c("A", "B", "A", "C", "B", "C", "A", "B", "C"))
})

# =============== test transform_digits =======================================
test_that("transform_digits works correcty", {
  # gets all negatives correctly
  expect_equal(transform_digits("0_0_0_0", c("A", "B", "C", "D")), "All negative")

  # gets single positives correctly
  expect_equal(transform_digits("0_0_1_0", c("A", "B", "C", "D")), "C")

  # gets double positives correctly
  expect_equal(transform_digits("1_0_1_0", c("A", "B", "C", "D")), "AC")

  # gets all positives correctly
  expect_equal(transform_digits("1_1_1_1", c("A", "B", "C", "D")), "ABCD")

  # check for differing lengths
  expect_error(
    transform_digits("1_1_1_1", c("A", "B", "C", "D", "E")),
    "Length of binary digits and input genes differ."
  )

  expect_error(
    transform_digits("1_1_1_1_1", c("A", "B", "C", "D")),
    "Length of binary digits and input genes differ."
  )

  # check for null input
  expect_error(
    transform_digits(NULL, c("A", "B", "C", "D")),
    "No binary digits given."
  )

  expect_error(
    transform_digits("1_1_1_1", NULL),
    "Length of binary digits and input genes differ."
  )
})

# =============== test extract_target_positive_colnames =======================
test_that("extract_target_positive_colnames correctly extracts relevant columns", {
  # Mock column names
  colnames <- c(
    "Env+", "Psi+", "Gag+", "Pol+",
    "Env+Psi+", "Psi+Pol+", "Gag+Pol+",
    "Concentration Env.Psi positive for target (copies/ul)",
    "Intact concentration Env.Psi (copies/ul)",
    "random_column"
  )

  # Test with exact match
  expect_equal(
    extract_target_positive_colnames("Env", colnames),
    c("Env+", "Env+Psi+")
  )

  expect_equal(
    extract_target_positive_colnames("Psi", colnames),
    c("Psi+", "Env+Psi+", "Psi+Pol+")
  )

  expect_equal(
    extract_target_positive_colnames("Pol", colnames),
    c("Pol+", "Psi+Pol+", "Gag+Pol+")
  )

  # Test with no match
  expect_equal(extract_target_positive_colnames("XYZ", colnames), character(0))

  # Test to ensure it excludes columns with spaces
  expect_false(any(grepl(" ", extract_target_positive_colnames("Env", colnames))))
  expect_false(any(grepl(" ", extract_target_positive_colnames("Psi", colnames))))
  expect_false(any(grepl(" ", extract_target_positive_colnames("Pol", colnames))))
})

# =============== test match_channel_gene =====================================
test_that("match_channel_gene correctly maps channels to genes", {
  # Mock input data
  df <- data.frame(
    `DyeName(s)` = c("FAM", "VIC", "Cy5", "ROX", "ATTO590"),
    Target = c("Env", "Psi", "Gag", "Pol", "RU5"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ch_dye <- c("Ch1" = "FAM", "Ch2" = "VIC", "Ch3" = "Cy5", "Ch5" = "ROX", "Ch6" = "ATTO590")

  # Test for 5 channels
  expect_equal(match_channel_gene(df, ch_dye, 5), c("Env+", "Psi+", "Gag+", "Pol+", "RU5+"))

  # Test for 4 channels (should exclude 5th dye)
  expect_equal(match_channel_gene(df, ch_dye, 4), c("Env+", "Psi+", "Gag+", "RU5+"))

  # Test invalid number of channels
  expect_error(match_channel_gene(df, ch_dye, 3), "Works only for 4 or 5 channels so far, but can be expanded easily.")

  # Test handling of NULL input
  expect_error(match_channel_gene(NULL, ch_dye, 5), "input df is not in right format")

  # Test missing columns in df
  df_invalid <- data.frame(Dye = c("FAM", "VIC"), stringsAsFactors = FALSE)
  expect_error(match_channel_gene(df_invalid, ch_dye, 5), "input df is not in right format")

  # Test multiple genes per dye (should trigger error)
  df_conflict <- df
  df_conflict$`DyeName(s)`[1] <- "VIC" # Simulate conflicting gene assignment
  expect_error(match_channel_gene(df_conflict, ch_dye, 5), "More than one gene per dye detected.")
})

# =============== test sum_target_positive_values ==============================
test_that("sum_target_positive_values correctly computes target sum", {
  mock_df <- data.frame(
    Target = c("Env", "Psi", "Gag"),
    `Env+` = c(5, NA, 3),
    `Psi+` = c(2, 6, 1),
    `Env+Psi+` = c(1, 2, 3),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_equal(sum_target_positive_values(mock_df[1, ]), 6)
  expect_equal(sum_target_positive_values(mock_df[2, ]), 8)
  expect_equal(sum_target_positive_values(mock_df[3, ]), 0)

  # Test missing Target column
  expect_error(sum_target_positive_values(mock_df[, -1]), "Data frame must contain a 'Target' column.")

  # Test empty data frame
  expect_equal(sum_target_positive_values(data.frame(Target = character())), 0)

  # Test no matching columns
  expect_equal(sum_target_positive_values(data.frame(Target = "XYZ")), 0)
})

# =============== test compute_groupwise_mean ==================================
test_that("compute_groupwise_mean correctly computes group means", {
  df <- data.frame(
    Well = rep(c("A1", "B1", "C1"), each = 2),
    Target = rep(c("Env", "Psi"), 3),
    Measurement = c(2, 3, NA, 6, 3, 9)
  )

  df_result <- compute_groupwise_mean(df, "Target", "Mean Measurement", "Measurement")

  expect_equal(df_result[["Mean Measurement"]], rep(c(2.5, 6), 3))
})

test_that("compute_groupwise_mean handles edge cases correctly", {
  df_empty <- data.frame(Well = character(), Target = character(), Measurement = numeric())
  expect_error(compute_groupwise_mean(df_empty, "Well", "Mean Measurement", "Old_Measurement"), "Specified columns must exist in the data frame.")

  df_non_numeric <- data.frame(Well = c("A1", "B1"), Target = c("Env", "Psi"), Measurement = c("high", "low"))
  expect_error(compute_groupwise_mean(df_non_numeric, "Well", "Mean Measurement", "Measurement"), "The column from which mean is computed must be numeric.")
})

# =============== test create_table ===========================================
# Sample input data
conf_mat <- data.frame(
  Well = c("A1", "B1", "C1", "D1"),
  `0_0_0_0` = c(10, 15, 20, 23),
  `1_0_0_0` = c(12, 14, 20, 19),
  `0_1_0_0` = c(21, 30, 12, 18),
  `0_0_1_0` = c(10, 11, 20, 13),
  `0_0_0_1` = c(2, 0, 2, 1),
  `1_1_0_0` = c(5, 5, 5, 1),
  `0_1_1_0` = c(3, 3, 5, 1),
  `1_0_1_0` = c(2, 3, 2, 1),
  `1_1_1_0` = c(2, 3, 1, 1),
  check.names = FALSE
)

dtQC <- data.frame(
  Well = rep(c("A1", "B1", "C1", "D1"), 4),
  Target = c(rep("Env", 4), rep("Psi", 4), rep("Gag", 4), rep("Pol", 4)),
  `DyeName(s)` = c(rep("VIC", 4), rep("Cy5", 4), rep("FAM", 4), rep("ATTO590", 4)),
  `Sample description 1` = rep("Sample1", 16),
  `Conc(copies/uL)` = rep(1, 16),
  `Accepted Droplets` = rep(c(8000, 5000, 9000, 7000), 4), # does not have to be the same, but simplicity
  Positives = c(10, 0, 2, 4, 5, 2, 6, 4, 1, 4, 5, 6, 3, 4, 2, 1),
  Negatives = rep(c(7000, 5000, 8800, 4000), 4), # does not have to be the same, but simplicity
  check.names = FALSE,
  stringsAsFactors = FALSE
)


tab1 <- data.frame(
  `Sample description 1` = c("Sample1", "Sample2"),
  `Mean concentration RPP30 + RPP30Shear (copies/uL)` = c(10, 10),
  `Mean unsheared` = c(0.5, 0.5),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

batch <- c("A1", "B1", "C1", "D1")
num_target <- 4
ch_dye <- list("Ch1" = "FAM", "Ch2" = "VIC", "Ch3" = "Cy5", "Ch5" = "ROX", "Ch6" = "ATTO590")
multi_pos <- get_multipos(genes = c("Env", "Psi", "Gag", "Pol"))
thresh <- 5000
tar_mio_factor <- c("Sample1" = 1.5, "Sample2" = 2.0)


# Test cases
test_that("create_table correctly processes batch data", {
  expect_warning(result <- create_table(
    c("A1", "B1", "C1", "D1"),
    conf_mat, dtQC, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1
  ))
  expect_type(result, "list")
  # Ensure added columns: 7 from dtQC, conf_mat, 8 added manually, and 9 per multiple positive
  expect_equal(ncol(result[[1]]), 7 + ncol(conf_mat) + 2 + 2 + 3 + 9 * length(multi_pos))
  expect_true("Target/Mio cells" %in% colnames(result[[1]]))
})

test_that("create_table handles H2O wells correctly", {
  dtQC_h2o <- dtQC
  dtQC_h2o$`Sample description 1` <- "H2O"
  result <- create_table(c("A1", "B1", "C1", "D1"), conf_mat, dtQC_h2o, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1)
  expect_equal(result[[2]], "h2o")
})

test_that("create_table filters wells below threshold", {
  expect_warning(result <- create_table(
    c("A1", "B1", "C1", "D1"), conf_mat, dtQC, num_target, ch_dye, multi_pos,
    10000, tar_mio_factor, tab1
  ))
  expect_equal(result[[2]], "empty")
})

test_that("create_table handles missing sample descriptions", {
  dtQC_missing <- dtQC
  dtQC_missing$`Sample description 1`[1] <- "UnknownSample"
  expect_error(
    create_table(c("A1", "B1", "C1", "D1"), conf_mat, dtQC_missing, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1),
    "Sample description not exactly matched"
  )
})

# =============== test create_tables ===========================================
# Input data
# Sample grouped data
grouped_data <- tibble(data.frame(
  `Sample description 1` = c("Sample1", "Sample2"),
  `Dilution Factor` = c(1, 1),
  RPP30 = c(NA, NA),
  shared_Wells = c("A1, B1", "C1, D1"),
  check.names = FALSE
))

# Sample input CSV data (mimicking preprocessed raw data)
in_csv <- data.frame(Well = c(rep("A1", 4), rep("B1", 4), rep("C1", 4), rep("D1", 4)),
                     Target1 = rep(c(1, 0, 0, 0), 4),
                     Target2 = rep(c(0, 1, 0, 0), 4),
                     Target3 = rep(c(0, 0, 1, 0), 4),
                     Target4 = rep(c(0, 0, 0, 1), 4),
                     Count = rep(10, 16))

# Sample dtQC data
dtQC <- data.frame(
  Well = rep(c("A1", "B1", "C1", "D1"), 4),
  Target = c(rep("Env", 4), rep("Psi", 4), rep("Gag", 4), rep("Pol", 4)),
  `DyeName(s)` = c(rep("VIC", 4), rep("Cy5", 4), rep("FAM", 4), rep("ATTO590", 4)),
  `Sample description 1` = c(rep("Sample1", 8), rep("Sample2", 8)),
  `Conc(copies/uL)` = rep(1, 16),
  `Accepted Droplets` = rep(c(8000, 5000, 9000, 7000), 4), # does not have to be the same, but simplicity
  Positives = c(10, 0, 2, 4, 5, 2, 6, 4, 1, 4, 5, 6, 3, 4, 2, 1),
  Negatives = rep(c(7000, 5000, 8800, 4000), 4), # does not have to be the same, but simplicity
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Sample additional parameters
tab1 <- data.frame(
  `Sample description 1` = c("Sample1", "Sample2"),
  `Mean concentration RPP30 + RPP30Shear (copies/uL)` = c(10, 10),
  `Mean unsheared` = c(0.5, 0.5),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

batch <- c("A1", "B1", "C1", "D1")
num_target <- 4
ch_dye <- list("Ch1" = "FAM", "Ch2" = "VIC", "Ch3" = "Cy5", "Ch5" = "ROX", "Ch6" = "ATTO590")
multi_pos <- get_multipos(genes = c("Env", "Psi", "Gag", "Pol"))
thresh <- 5000
tar_mio_factor <- c("Sample1" = 1.5, "Sample2" = 2.0)

# test
test_that("create_tables processes valid input correctly", {
  # Mock `create_table` to return controlled values
  mock_create_table <- mock(
    list(data.frame(Well = c("A1", "B1"), `Processed Value` = c(5, 6)), "tab"),
    list(data.frame(Well = c("C1", "D1"), `Processed Value` = c(3, 4)), "tab")
  )
  
  stub(create_tables, "create_table", mock_create_table)
  
  result <- create_tables(grouped_data, in_csv, dtQC, ch_dye, multi_pos, thresh, tar_mio_factor, tab1)
  
  # Extract results
  output_tables <- result[[1]]
  conf_mats <- result[[2]]
  h2o_tables <- result[[3]]
  
  # Validate output
  expect_length(output_tables, 2)  # Two groups should be processed
  expect_true(all(sapply(output_tables, is.data.frame)))
  
  expect_length(conf_mats, 2)  # Two confusion matrices should be created
  expect_true(all(sapply(conf_mats, is.data.frame)))
  
  expect_length(h2o_tables, 0)  # No water wells in this case
})

test_that("create_tables handles wells with <4 targets", {
  in_csv_3 <- data.frame(Well = c(rep("A1", 4), rep("B1", 4), rep("C1", 4), rep("D1", 4)),
                       Target1 = rep(c(1, 0, 0, 0), 4),
                       Target2 = rep(c(0, 1, 0, 0), 4),
                       Target3 = rep(c(0, 0, 1, 0), 4),
                       Count = rep(10, 16))

  
  expect_warning(
    create_tables(grouped_data, in_csv_3, dtQC, ch_dye, multi_pos, thresh, tar_mio_factor, tab1),
    "Less than 4 Targets detected. Skipping."
  )
})

test_that("create_tables correctly processes water control wells", {
  dtQC_2 <- data.frame(
    Well = c("A1", "B1"),
    `Sample description 1` = c("H2O", "H2O")  # Water wells
  )

  
  mock_create_table_2 <- mock(list(data.frame(Well = c("A1", "B1"), `Processed Value` = c(NA, NA)), "h2o"), cycle = TRUE)
  
  stub(create_tables, "create_table", mock_create_table_2)
  
  result <- create_tables(grouped_data, in_csv, dtQC_2, ch_dye, multi_pos, thresh, tar_mio_factor, tab1)
  
  # Expect water wells to be stored in h2o_tables
  expect_length(result[[1]], 0)  # No standard output tables
  expect_length(result[[2]], 2)  # One confusion matrix
  expect_length(result[[3]], 2)  # Water table should be stored
})

