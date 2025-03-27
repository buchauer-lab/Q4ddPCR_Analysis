# ======================== Beauty Output =======================================

#' Create analysis tables
#'
#' Create an output sheet that collects the most important information.
#' @param table Data frame computed in create_tables().
#' @param multi_pos The possible multiple positives.
#' @return A named vector as blueprint for a xlsx sheet.
get_output_sheet <- function(table, multi_pos, grouped_data) {
  # create output of single positives
  sample <- unique(table$`Sample description 1`)
  output_sheet <- c(
    "Title" = paste(unique(table$Well), collapse = ", "),
    "Sample" = sample,
    "Number of cells analysed" = pull(unique(round(grouped_data[grouped_data$`Sample description 1` == sample, "Mean cells per reaction"]))),
    "Shearing Index" = pull(unique(round(grouped_data[grouped_data$`Sample description 1` == sample, "Shearing index"], 2))),
    "Total HIV-DNA (4-based) (copies/Mio cells)" = round(unique(table[["total HIV DNA/Mio cells"]])),
    "Total HIV-DNA (Env.Psi) (copies/Mio cells)" = round(unique(table[["total HIV DNA/Mio cells (Env.Psi)"]])),
    "Gag/Mio cells" = round(unique(table[grepl("Gag", table$Target), "Mean Target/Mio cells"])),
    "Pol/Mio cells" = round(unique(table[grepl("Pol", table$Target), "Mean Target/Mio cells"])),
    "Psi/Mio cells" = round(unique(table[grepl("Psi", table$Target), "Mean Target/Mio cells"])),
    "Env/Mio cells" = round(unique(table[grepl("Env", table$Target), "Mean Target/Mio cells"])),
    "Psi-defective" = round(unique(table[grepl("Psi", table$Target), "Mean Target/Mio cells"]) -
      unique(na.omit(table[["intact provirus/Mio cells Psi.Env, corrected for shearing"]]))),
    "Env-defective" = round(unique(table[grepl("Env", table$Target), "Mean Target/Mio cells"]) -
      unique(na.omit(table[["intact provirus/Mio cells Psi.Env, corrected for shearing"]])))
  )

  # loop over multiple positives and add to output
  for (multip in multi_pos) {
    name0 <- paste0("Intact concentration ", paste(multip, collapse = "."), " (copies/ul)")
    name1 <- paste0("intact provirus/Mio cells ", paste(multip, collapse = "."))
    name2 <- paste0(name1, ", corrected for shearing")

    x0 <- round(unique(na.omit(table[[name0]])), 2)
    x1 <- round(unique(na.omit(table[[name1]])))
    x2 <- round(unique(na.omit(table[[name2]])))

    output_sheet[[name0]] <- ifelse(length(x0) != 0,
      x0,
      NA
    )
    output_sheet[[name1]] <- ifelse(length(x1) != 0,
      x1,
      NA
    )
    output_sheet[[name2]] <- ifelse(length(x2) != 0,
      x2,
      NA
    )
  }
  return(output_sheet)
}

# =================== Save to xlsx file =======================================

#' Write data to xlsx
#'
#' Write finished analysis to a xlsx file
#' @param output_tables List of data frames. Output from create_tables().
#' @param conf_mats The computed confusion matrices.
#' @param tab1 Data frame created in create_household_table().
#' @param output_file File name (path) of the output file.
#' @param h2o_table Data frame with H2O samples created in create_tables().
#' @param multi_pos The possible multiple positives.
#' @return No return value. Creates an xlsx file as output_file.
#' @export
write_output_file <- function(output_tables, conf_mats, tab1, output_file, h2o_table, multi_pos) {
  l <- lapply(output_tables, get_output_sheet, multi_pos, grouped_data)
  output_sheet <- do.call(rbind.data.frame, l)

  output_list <- append(list(output_sheet), output_tables)
  output_list <- append(append(output_list, list(tab1)), h2o_table)

  openxlsx::write.xlsx(output_list, output_file)
}
