#' Merging information of both input files
#'
#' This function creates one table based on the two input files: in_xlsx (as dtQC
#' or data_table in tutorial) and in_csv (as confusion matrix in tutorial).
#' Merging is based on wells, information from shearing table will be added.
#' @param data_table Information coming from in_xlsx file, wells with targets only
#' @param conf_mat Output from create_confusion_matrix
#' @param shear_table Output from compute_shearing_factor
#' @import tibble
#' @return Combined table, conserving all information
#' @export
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

#' Compute target mean
#'
#' This function computes the number of targets (genes) detected per million cells
#' and the mean thereof within the wells of one group
#' @param tab Output from merge_tables
#' @import dplyr, tidyr
#' @return Updated table
#' @export
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
#' Compute the possible combinations of targets (multiplets) based on the targets
#' in the experiment.
#' @param targets List of targets used in the experiment.
#' @return List of possible combinations of input genes
#' @export
get_multipos <- function(targets) {
  multi_pos <- list()

  for (i in 2:length(targets)) {
    # Get the combinations of length i
    comb <- combn(targets, i, simplify = FALSE)
    # Append the combinations to the list
    multi_pos <- c(multi_pos, comb)
  }

  return(multi_pos)
}



#' Create a confusion matrix
#'
#' Create a confusion matrix showing which targets and multiplets (combination of
#' targets) are found in which wells
#' @param df in_csv subsetted for wells containing targets
#' @param data_table dtQC subsetted for wells containing target
#' @param ch_dye named list with used channels as names and respective dyes as values
#' @param target_channel named list with Target (in_cs columns) as names and channels as values
#' @import dplyr, tidyr, tibble
#' @return Confusion matrix
#' @export
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
