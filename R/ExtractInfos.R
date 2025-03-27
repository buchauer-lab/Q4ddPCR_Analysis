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

#' Convert digits into genes
#'
#' Convert the binary digits of into the genes they represent (i.e., match
#' target from csv file to genes).
#' @param combination String with binaries (0 and 1) separated by underscores
#' ("_").
#' @param genes Vector of genes that were targeted for.
#' @return Vector containing the transformed strings.
transform_digits <- function(combination, genes) {
  # check for input
  if (is.null(combination)) {
    stop("No binary digits given.")
  }

  # Split combination into individual digits
  digits <- strsplit(combination, "_")[[1]]

  # check for length mismatch
  if (length(digits) != length(genes)) {
    stop("Length of binary digits and input genes differ.")
  }
  # Initialize empty vector to store object names
  positive_objects <- character(0)

  # Iterate over each digit and add corresponding object name if positive
  for (i in 1:length(digits)) {
    if (as.integer(digits[i]) == 1) {
      positive_objects <- c(positive_objects, genes[i])
    }
  }

  # Concatenate object names
  transformed_combination <- paste(positive_objects, collapse = "")
  if (transformed_combination == "") {
    return("All negative")
  }
  return(transformed_combination)
}

#' Extract target positive column names
#'
#' Extract those column names that are positive for a specified target (gene).
#' @param target The target the column names should be positive for.
#' @param colnames The column names of a data frame.
#' @return Vector containing the gene positive column names.
extract_target_positive_colnames <- function(target, colnames) {
  cols <- grep(target, colnames, value = T)
  # remove spaces to not include column names that are added but where not part of the original data
  return(grep(" ", cols, value = T, invert = T))
}

#' Match measurement channel to gene
#'
#' Match the Channels that were measured with the genes using the dye information.
#' @param df Data frame with columns "DyeName(s)" and "Target".
#' @param ch_dye Named vector matching channels (names) with dyes (values).
#' @param num_ch The number of channels used in total (must be 4 or 5).
#' @import dplyr
#' @return Vector of genes.
match_channel_gene <- function(df, ch_dye, num_ch) {
  # check input format
  if (is.null(df) || !all(c("DyeName(s)", "Target") %in% colnames(df))) {
    stop("input df is not in right format")
  }
  # loop over dyes and match to gene
  genes <- c()
  for (dye in ch_dye) {
    gene <- df[df$"DyeName(s)" == dye, "Target"]
    if (length(unique(gene)) > 1) {
      stop("More than one gene per dye detected.")
    }
    gene_tmp <- substr(gene[1], nchar(gene[1]) - 2, nchar(gene[1]))
    genes <- append(genes, paste0(gene_tmp, "+"))
  }
  if (num_ch == 4) {
    return(genes[c(1, 2, 3, 5)])
  } else if (num_ch == 5) {
    return(genes)
  } else {
    stop("Works only for 4 or 5 channels so far, but can be expanded easily.")
  }
}

#' Compute total target positive values
#'
#' Compute sum of all target positive values (including doublets, triplets, etc.)
#' @param df Data frame with experiment information.
#' @import dplyr
#' @return Numeric: sum of all target positive values.
sum_target_positive_values <- function(df) {
  # check for column
  if (is.null(df) || !"Target" %in% names(df)) {
    stop("Data frame must contain a 'Target' column.")
  }
  # type check
  if (!is.character(df[["Target"]])) {
    stop("Target column must be of type character.")
  }
  # get Target
  target <- substr(df["Target"], (nchar(df["Target"]) - 2), nchar(df["Target"]))
  cols <- extract_target_positive_colnames(target, names(df))
  return(sum(as.numeric(df[cols]), na.rm = T))
}

#' Compute group-wise mean
#'
#' Compute the mean of a specified group of Wells for a certain column and add
#' as new column
#' @param df Data frame with experiment information.
#' @param group_by Column(s) by which the data should be grouped.
#' @param new_col The name of the new column.
#' @param old_col The name of the column the mean should be computed from.
#' @import dplyr
#' @return df updated by new_col.
compute_groupwise_mean <- function(df, group_by, new_col, old_col) {
  # check input format
  if (!all(c(old_col, group_by) %in% colnames(df))) {
    stop("Specified columns must exist in the data frame.")
  }
  
  if (!is.numeric(df[[old_col]])) {
    stop("The column from which mean is computed must be numeric.")
  }
  # compute mean (by group)
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

# TODO: Can probably be deleted, but I am to afraid yet
## ' Get relevant subsets
#'
#' Compute the possible subsets of a multiple positive target.
# #' @param multi_pos The multiple positive target.
# #' @return List of the possible subsets.
# #' @import dplyr
#get_subset <- function(multi_pos) {
#  # Initialize a list to store the subset relationships
#  subset_list <- vector("list", length(multi_pos))
#  names(subset_list) <- seq_along(multi_pos)

#  # Check which combinations are subsets of which other combinations
#  for (i in 1:(length(multi_pos) - 1)) {
#    subset_list[[i]] <- list()
#    for (j in 1:(length(multi_pos) - 1)) {
#      if (i != j) {
#        if (all(multi_pos[[i]] %in% multi_pos[[j]])) {
#          subset_list[[i]] <- c(subset_list[[i]], j)
#        }
#      }
#    }
#  }
#  subset_list <- lapply(subset_list, unlist)
#  return(subset_list)
#}

#' Create analysis table for batch
#'
#' Combine input data and compute analysis parameters for a group of wells (batch).
#' @param batch The wells which are analysed together.
#' @param conf_mat Confusion matrix (well vs. positives) for the batch.
#' @param dtQC dtQC table created in the Readers/read_files function.
#' @param tab1 Table with RPP30(Shear) information.
#' @param num_target The total number of targets.
#' @param ch_dye Named list to match channel (names) with dyes (values).
#' @param multi_pos The possible multiple positives.
#' @param thresh The minimum number of accepted droplets in well to be included
#' in the analysis.
#' @param tar_mio_factor Factor to multiply Concentration with to obtain
#' Target/Mio cells
#' @param tab1 Data frame created in create_household_table().
#' @import dplyr
#' @return List with two elements: a data frame and a sting indicating the nature
#' of the data sample (empty, H2O control, data)
create_table <- function(batch, conf_mat, dtQC, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1) {
  # create table of respective wells
  tab <- merge(dtQC[dtQC$Well %in% batch, 1:8], conf_mat[conf_mat$Well %in% batch, ], by = "Well", all = T)
  tab <- Filter(function(x) !all(is.na(x)), tab)
  
  if(length(unique(tab$`Sample description 1`)) > 1){
    stop("Sample description not exactly matched")
  }

  # get correct tar_mio_factor
  tar_mio <- tar_mio_factor[unique(tab$`Sample description 1`)]

  # change names (from 1_0_0_0 to Gag+ etc)
  markers <- match_channel_gene(tab, ch_dye, num_target)
  names(tab)[grepl("_", names(tab))] <- unname(sapply(grep("_", names(tab), value = T), transform_digits, markers))

  # remove wells with less than x (7500) accepted droplets
  tab <- sufficient_droplets(tab, thresh)

  # check if any Well left after thresholding droplets
  if (nrow(tab) == 0) {
    return(list("", "empty"))
  }

  # do not compute values for H2O Wells
  if (unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")) {
    return(list(tab, "h2o"))
  }

  # add mean concentration and mean unsheared

  # check if sample description is in mean concentration household array
  if (unique(tab$`Sample description 1`) %in% tab1$`Sample description 1`) {
    tab$`Mean concentration RPP30 + RPP30Shear (copies/uL)` <- unique(tab1[tab1$`Sample description 1` == unique(tab$`Sample description 1`), "Mean concentration RPP30 + RPP30Shear (copies/uL)"])[[1]]
    # mean_conc_household[unique(tab$`Sample description 1`)]
    tab$`Mean unsheared` <- unique(tab1[tab1$`Sample description 1` == unique(tab$`Sample description 1`), "Mean unsheared"])[[1]]
  } else {
    # if not, stop procedure
    stop("Sample description not exactly matched in names for Mean RPP30(Shear) concentration and Mean unsheared.
         Make sure to match Sample descriptions for RPP samples and HIV-target samples.")
  }

  # compute Target per million cells
  tab <- tab %>%
    mutate(`Target/Mio cells` = tar_mio * 10^6 * `Conc(copies/uL)` / (`Mean concentration RPP30 + RPP30Shear (copies/uL)`))

  # compute mean target per mio cells for each "group" (e.g. A09, B09, C09, D09)
  tab <- compute_groupwise_mean(
    tab, c("Target"), "Mean Target/Mio cells",
    "Target/Mio cells"
  )

  # compute values for multiple positives for specified combinations
  for (multip in multi_pos) {
    tab <- get_multi_pos(tab, multip, tar_mio)
  }

  tab <- compute_total_HIV(tab)

  tab <- compute_total_HIV_envPsi(tab)

  return(list(tab, "tab"))
}



# ====================== Execute ===============================================
#' Create analysis tables
#'
#' Combine input data and do analysis.
#' @param grouped_data Data frame that gives information per group of wells.
#' Computed in define_groups()
#' @param in_csv Data frame with information from csv input file. Output from
#' read_csv().
#' @param dtQC dtQC table created in the Readers/read_files function.
#' @param ch_dye Named list to match channel (names) with dyes (values).
#' @param multi_pos The possible multiple positives.
#' @param thresh The minimum number of accepted droplets in well to be included
#' in the analysis.
#' @param tar_mio_factor Factor to multiply Concentration with to obtain
#' Target/Mio cells
#' @param tab1 Data frame created in create_household_table().
#' @return 3 Lists containing the (1) data analysis, (2) confusion matrices, (3)
#' the H2O samples.
#' @export
create_tables <- function(grouped_data, in_csv, dtQC, ch_dye, multi_pos, thresh, tar_mio_factor, tab1) {
  # initialize empty lists to save output
  output_tables <- list()
  h2o_tables <- list()
  conf_mats <- list()

  # loop over groups and run functions
  for (group in apply(
    grouped_data[is.na(grouped_data$RPP30), "shared_Wells"],
    1,
    function(x) (unname(strsplit(x, ", ")))
  )) {
    group <- unlist(group)
    # get relevant section of csv
    sub_in_csv <- tibble(in_csv[in_csv$Well %in% group, ])
    # remove columns with all na
    sub_in_csv <- sub_in_csv[, !apply(is.na(sub_in_csv), 2, all)]
    # compute confusion matrix
    # TODO: before conf_mat was computed once for all wells with the same targets
    # TODO: append appropriate conf matrices at the end
    conf_mat <- tidyr::pivot_wider(sub_in_csv,
      id_cols = Well,
      names_from = grep("Target", names(sub_in_csv), value = T),
      values_from = Count,
      values_fill = 0
    )
    # get number of targets
    num_target <- length(grep("Target", names(sub_in_csv), value = T))
    # check number of targets
    if (num_target < 4) {
      warning("Less than 4 Targets detected. Skipping. If this is expexted
              behaviour, please contact Mark :)")
      next
    }

    # compute output tables (equivalent to analysis sheet in manual analysis)

    # create output table
    out_table <- create_table(group, conf_mat, dtQC, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1)

    # save to water table if water control
    if (out_table[[2]] == "h2o") {
      h2o_tables[[length(h2o_tables) + 1]] <- out_table[[1]]
    } else if (out_table[[2]] == "empty") {
      # skip if empty (everything removed during QC)
      next
    } else {
      output_tables[[length(output_tables) + 1]] <- out_table[[1]]
    }


    # save confusion matrix
    conf_mats[[length(conf_mats) + 1]] <- conf_mat
  }
  return(list(output_tables, conf_mats, h2o_tables))
}
