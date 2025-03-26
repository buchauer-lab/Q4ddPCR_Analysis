library(testthat)
source("../R/MultipPositives.R")

#############
test_that("get_multi_pos works correctly", {
  set.seed(0)

  # Mock input data
  df <- data.frame(
    Well = rep(c("A1", "B1", "C1", "D1"), 4),
    Target = c(rep("Env", 4), rep("Psi", 4), rep("Gag", 4), rep("Pol", 4)),
    `Conc(copies/µL)` = rep(1, 16),
    `Env+` = rep(c(6, 7, 4, 4), 4),
    `Psi+` = rep(c(8, 8, 1, 3), 4),
    `Pol+` = rep(c(1, 2, 3, 4), 4),
    `Gag+` = rep(c(1, 4, 3, 2), 4),
    `Env+Psi+` = rep(c(5, 7, 1, 2), 4),
    `Psi+Pol+` = rep(c(0, 1, 1, 1), 4),
    `Gag+Pol+` = rep(c(0, 1, 1, 1), 4),
    `Env+Psi+Pol+` = rep(c(0, 1, 0, 1), 4),
    `Mean concentration RPP30 + RPP30Shear (copies/µL)` = rep(100, 16),
    `Mean unsheared` = rep(0.5, 16),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tar_mio_factor <- 1
  genes <- c("Env", "Psi")

  # Mock functions
  sum_target_positive_values <- function(df) {
    target <- substr(df["Target"], (nchar(df["Target"]) - 2), nchar(df["Target"]))
    cols <- extract_target_positive_colnames(target, names(df))
    return(sum(as.numeric(df[cols]), na.rm = T))
  }

  compute_groupwise_mean <- function(df, group_by, new_col, old_col) {
    if (grepl("Mean|Intact", new_col)) {
      df <- df %>%
        group_by_at(vars(!!!group_by)) %>%
        mutate(!!sym(new_col) :=
          mean(!!sym(old_col), na.rm = TRUE)) %>%
        mutate(!!sym(paste0("STD ", new_col)) :=
          sd(!!sym(old_col), na.rm = TRUE)) %>%
        ungroup()
    } else {
      df <- df %>%
        group_by_at(vars(!!!group_by)) %>%
        mutate(!!sym(new_col) :=
          mean(!!sym(old_col), na.rm = TRUE)) %>%
        ungroup()
    }

    return(df)
  }

  # Run function
  df_result <- get_multi_pos(df, genes, tar_mio_factor)

  # Check that expected columns are created
  expect_true("Concentration Env.Psi positive for target (copies/µl)" %in% names(df_result))
  expect_true("Intact concentration Env.Psi (copies/µl)" %in% names(df_result))
  expect_true("intact provirus/Mio cells Env.Psi" %in% names(df_result))
  expect_true("intact provirus/Mio cells Env.Psi, corrected for shearing" %in% names(df_result))

  # Check that results are numeric
  expect_type(df_result$`Concentration Env.Psi positive for target (copies/µl)`, "double")
  expect_type(df_result$`Intact concentration Env.Psi (copies/µl)`, "double")
  expect_type(df_result$`intact provirus/Mio cells Env.Psi`, "double")
  expect_type(df_result$`intact provirus/Mio cells Env.Psi, corrected for shearing`, "double")

  # Check na values are correctly computed
  expect_true(sum(is.na(df_result$`Concentration Env.Psi positive for target (copies/µl)`)) == 8)
  expect_true(sum(is.na(df_result$`Intact concentration Env.Psi (copies/µl)`)) == 8)
  expect_true(sum(is.na(df_result$`intact provirus/Mio cells Env.Psi`)) == 8)
  
  # get warning if not all genes appear as multiplets
  expect_warning(get_multi_pos(df, c("Env", "Gag"), tar_mio_factor),
                 "No multiple positives were found in data for:EnvGag",
                 fixed=TRUE)

  # Check values are correct
  expect_equal(unique(na.omit(df_result$`Intact concentration Env.Psi (copies/µl)`)),
    0.4041948,
    tolerance = 1e-06
  )
  expect_equal(unique(na.omit(df_result$`intact provirus/Mio cells Env.Psi, corrected for shearing`)),
    8083.896,
    tolerance = 1e-06
  )
})
