# this is the beginning of the create_table function (maybe can be done for all at once)
merge_tables <- function(data_table, conf_mat, shear_table){
  tab <- merge(data_table, conf_mat, by = "Well", all = T)
  tab <- Filter(function(x) !all(is.na(x)), tab)
  
  # add mean RPP concentration and mean unsheared to table
  if (all(unique(tab$`Sample description 1`) %in% shear_table$`Sample description 1`)) {
    mean_rpp_conc <- tibble::deframe(unique(shear_table[,c("Sample description 1", "Mean concentration RPP30 + RPP30Shear (copies/uL)")]))
    mean_unsheared <- tibble::deframe(unique(shear_table[,c("Sample description 1", "Mean unsheared")]))
    tab$`Mean concentration RPP30 + RPP30Shear (copies/uL)` <- mean_rpp_conc[tab$`Sample description 1`]
    tab$`Mean unsheared` <- mean_unsheared[tab$`Sample description 1`]
  } else {
    # if not, stop procedure
    stop("revise this")
  }
  return(tab)
}

compute_target_means <- function(tab){
  # compute target per million cells
  tab <- tab %>%
    mutate(`Target/Mio cells` = tar_mio_factor[`Sample description 1`] * 10^6 * `Conc(copies/uL)` / (`Mean concentration RPP30 + RPP30Shear (copies/uL)`))
  
  # compute mean target per million cells (mean per group and per target)
  tab <- tab %>%
    group_by(group_id, Target) %>%
    mutate(`Mean Target/Mio cells` = mean(`Target/Mio cells`, na.rm = TRUE)) %>%
    mutate(`SD Target/Mio cells` = sd(`Target/Mio cells`, na.rm = TRUE)) %>%
    ungroup()
  
  return(tab)
}

# ====================== Define functions to calculate table ==================

#' Get multiple positive combinations
#'
#' Compute the possible combinations of genes based on the genes in an experiment.
#' @param genes List of genes used in the experiment.
#' @return List of possible combinations
#' @export
get_multipos <- function(genes) {
  multi_pos <- list()

  for (i in 2:length(genes)) {
    # Get the combinations of length i
    comb <- combn(genes, i, simplify = FALSE)
    # Append the combinations to the list
    multi_pos <- c(multi_pos, comb)
  }

  return(multi_pos)
}




create_confusion_matrix <- function(df, data_table, ch_dye, target_channel){
  # df is in_csv subsetted for data wells
  
  # map correct gene names
  ## rename Targets to corresponding gene names
  target_dye_map <- setNames(ch_dye[target_channel], names(target_channel))
  dye_gene_map <- tibble::deframe(unique(data_table[,c("DyeName(s)", "Target")]))
  ## TODO: This only works when gene is at end of name (e.g., SilicanoEnv)
  dye_gene_map <- unlist(lapply(dye_gene_map, function(x){paste0(substr(x, nchar(x)-2, nchar(x)),'+')}))
  ## get index of csv Target columns
  idx <- grep("Target", names(df))
  ## rename Target columns
  names(df)[idx] <- dye_gene_map[target_dye_map[names(df)[idx]]]
  
  # add the names as column (e.g., 'Psi+' to 1_0_0_0)
  # TODO: should not be based on indices
  df_long <- df[,1:6] %>%
    rowwise() %>%
    mutate(
      name = {
        cols <- names(across(2:5))[c_across(2:5) == 1]
        if (length(cols) == 0) "all_negative" else paste(cols, collapse = "")
      }
    ) %>%
    ungroup()
  # compute confusion table
  conf_mat <- df_long %>%
    tidyr::pivot_wider(
      id_cols = Well,
      names_from = name,
      values_from = Count,
      values_fill = 0
    )
  return(conf_mat)
}
