# ======================== Beauty Output =======================================

#' Create analysis tables
#'
#' Create an output sheet that collects the most important information.
#' @param table Data frame computed in create_tables().
#' @param multi_pos The possible multiple positives.
#' @param shear_table Output from compute_shearing_factor.
#' @return A named vector as blueprint for a xlsx sheet.
#' @keywords internal
get_output_sheet <- function(table, multi_pos, shear_table) {
  # create output of single positives
  sample <- unique(table$`Sample description 1`)
  output_sheet <- c(
    "Title" = paste(unique(table$Well), collapse = ", "),
    "Sample" = sample,
    "Number of cells analysed" = unique(round(shear_table[shear_table$`Sample description 1` == sample, "Mean cells per reaction"])),
    "Shearing Index" = pull(unique(round(shear_table[shear_table$`Sample description 1` == sample, "Shearing index"], 2))),
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
  
  # save multiplets
  out_names <- grep("^Intact concentration|^intact provirus/Mio cells", names(table), value=T)
  out_names1 <- grep("^Intact concentration", names(table), value=T)
  out_names2 <- grep("^intact provirus/Mio cells", names(table), value=T)
  out_tab <- unique(cbind(round(table[,out_names1],2),round(table[,out_names2]))[,out_names])
  out_list <- lapply(out_tab, function(col) unique(col[!is.na(col)]))
  output_sheet <- c(output_sheet, out_list)

  return(output_sheet)
}

# =================== Save to xlsx file =======================================

#' Write data to xlsx
#'
#' Write finished analysis to a xlsx file
#' @param output_tables List of data frames. Output from create_tables().
#' @param conf_mats The computed confusion matrices.
#' @param output_file File name (path) of the output file.
#' @param h2o_table Data frame with H2O samples created in create_tables().
#' @param multi_pos The possible multiple positives.
#' @param shear_table Output from compute_shearing_factor.
#' @return No return value. Creates an xlsx file as output_file.
#' @export
write_output_file <- function(output_tables, conf_mats, output_file, h2o_table, multi_pos, shear_table) {
  l <- lapply(output_tables, get_output_sheet, multi_pos, shear_table)
  output_sheet <- do.call(rbind.data.frame, l)

  output_list <- append(list(output_sheet), output_tables)
  output_list <- append(append(output_list, list(shear_table)), h2o_table)

  openxlsx::write.xlsx(output_list, output_file)
}
