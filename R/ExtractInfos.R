# branch dev

# extract relevant information from csv and xlsx file

# load necessary packages
library(readxl)
library(tidyr)
library(dplyr)
source("R/MultipPositives.R")
source("R/CreateHouseholdTable.R")

# ====================== Define functions to calculate table ==================

#' Get multiple positive combinations
#' 
#' Compute the possible combinations of genes based on the genes in an experiment.
#' @param genes List of genes used in the experiment.
#' @return List of possible combinations
#' @export
get_multipos <- function(genes){
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
#' @export
transform_digits <- function(combination, genes) {
  
  # Split combination into individual digits
  digits <- strsplit(combination, "_")[[1]]

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
  if(transformed_combination == ""){
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
#' @export
extract_target_positive_colnames <- function(target, colnames){
  cols <- grep(target, colnames, value=T)
  return(grep(" ", cols, value=T, invert = T))
}

#' Match measurement channel to gene
#' 
#' Match the Channels that were measured with the genes using the dye information.
#' @param df Data frame with columns "DyeName(s)" and "Target".
#' @param ch_dye Named vector matching channels (names) with dyes (values).
#' @param num_ch The number of channels used in total (must be 4 or 5).
#' @return Vector of genes.
#' @export
match_channel_gene <- function(df, ch_dye, num_ch){
  genes <- c()
  for (dye in ch_dye){
    gene <- df[df$"DyeName(s)" == dye, "Target"]
    if(length(unique(gene)) > 1){
      stop("More than one gene per dye detected.")
    }
    gene_tmp <- substr(gene[1], nchar(gene[1])-2, nchar(gene[1]))
    genes <- append(genes, paste0(gene_tmp, "+"))
  }
  if(num_ch == 4){
    return(genes[c(1,2,3,5)])
  }
  else if(num_ch == 5){
    return(genes)
  } else{
    stop("Works only for 4 or 5 channels so far, but can be expanded easily.")
  }
}

#' Compute total target positive values
#' 
#' Compute sum of all target positive values (including doublets, triplets, etc.)
#' @param df Data frame with experiment information.
#' @return Numeric: sum of all target positive values.
#' @export
sum_target_positive_values <- function(df){
  target <- substr(df["Target"], (nchar(df["Target"])-2), nchar(df["Target"]))
  cols <- extract_target_positive_colnames(target, names(df))
  return(sum(as.numeric(df[cols]), na.rm=T))
}

#' Compute group-wise mean
#' 
#' Compute the mean of a specified group of Wells for a certain column and add
#' as new column
#' @param df Data frame with experiment information.
#' @param group_by Column(s) by which the data should be grouped.
#' @param new_col The name of the new column.
#' @param old_col The name of the column the mean should be computed from.
#' @return df updated by new_col.
#' @export
compute_groupwise_mean<- function(df, group_by, new_col, old_col){
  if(grepl("Mean|Intact", new_col)){
    df <- df %>%
      group_by_at(vars(!!!group_by)) %>%
      mutate(!!sym(new_col) :=
               mean(!!sym(old_col), na.rm = TRUE)) %>%
      mutate(!!sym(paste0("STD ", new_col)) :=
               sd(!!sym(old_col), na.rm = TRUE)) %>%
      ungroup()
  } else{
    df <- df %>%
      group_by_at(vars(!!!group_by)) %>%
      mutate(!!sym(new_col) :=
               mean(!!sym(old_col), na.rm = TRUE)) %>%
      ungroup()
  }

  return(df)
}

#' Get relevant subsets
#' 
#' Compute the possible subsets of a multiple positive target.
#' @param multi_pos The multiple positive target.
#' @return List of the possible subsets.
#' @export
get_subset <- function(multi_pos){
  # Initialize a list to store the subset relationships
  subset_list <- vector("list", length(multi_pos))
  names(subset_list) <- seq_along(multi_pos) 
  
  # Check which combinations are subsets of which other combinations
  for (i in 1:(length(multi_pos)-1)) {
    subset_list[[i]] <- list()
    for (j in 1:(length(multi_pos)-1)) {
      if (i != j) {
        if (all(multi_pos[[i]] %in% multi_pos[[j]])) {
          subset_list[[i]] <- c(subset_list[[i]], j)
        }
      }
    }
  }
  subset_list <- lapply(subset_list, unlist)
  return(subset_list)
}

#' Compute total HIV content
#' 
#' Compute the total HIV content as proposed by Rachel Scheck.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @param multi_pos List of multiple positives. #TODO: check again.
#' @return Data frame df updated with new columns.
#' @export
compute_total_HIV <- function(df, multi_pos){
  # get length of multiple positives
  len_pos <- unlist(lapply(multi_pos, length))
  
  # find which multiple positives are a subset of another one
  subsets <- get_subset(multi_pos)

  # check that there are multiple positives that need to be corrected
  if(max(len_pos) < 2){
    warning("Maximum number of genes in multiple positives is 2. Total HIV DNA
            per Mio cells will not be computed")
    return(df)
  }
  
  # get shear corrected all positive value
  quad_value_shear_corrected <- df[[paste0("intact provirus/Mio cells ",
                                            paste(multi_pos[[which.max(len_pos)]], collapse = "."),
                                            ", corrected for shearing")]]
  # loop over the length of combinations
  for (i in (max(len_pos)-1):2){
    which_multipos <- which(len_pos == i)

   # loop over all multi positives
    for(comb in which_multipos){
      # get gene names
      genes <- multi_pos[[comb]]
      # get name for new variable
      name <- paste0("Higher multipos corrected ", paste(genes, collapse = "."))
      # get value that needs to be corrected
      original_value <- df[[paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))]]
      # correct by all positive value
      new_value = original_value - quad_value_shear_corrected
      # correct by subsets if necessary
      if(!is.null(subsets[[as.character(comb)]])){
        for(subset in subsets[[as.character(comb)]]){
          # substract subset value
          subset_value <- df[[paste0("Higher multipos corrected ", paste(multi_pos[[subset]], collapse = "."))]]
          #print(subset_value)
          new_value = new_value - subset_value
        }
      }
      # save new values to table

      df[[name]] <- unique(na.omit(new_value))
    }
    
  }
  # correct for single positives
  target <- substr(df$Target, nchar(df$Target)-2, nchar(df$Target))

  # compute for each column
  subtract_multipos <- lapply(unique(target), function(tar){
    cols.of.int <-  grep(tar, grep("Higher multipos corrected", names(df), value=T), value=T)
    res <- apply(df[grepl(tar, df$Target),cols.of.int], 1, sum, na.rm=TRUE)
    return(res)
  })

  df[["Multipos corrected Mean Target/Mio cells"]] <- df[["Mean Target/Mio cells"]] - quad_value_shear_corrected - unlist(subtract_multipos)
  
  df[["Higher multipos corrected sum for all targets/Mio cells"]] <- 4 * mean(df[["Multipos corrected Mean Target/Mio cells"]], na.rm=TRUE)
  
  # TODO
  df[["total HIV DNA/Mio cells"]] <-  quad_value_shear_corrected + apply(df[,grep("Higher multipos corrected", names(df), value=TRUE)], 1, sum, na.rm=T)
  return(df)
}

#' Compute total HIV content
#' 
#' Estimate the total HIV content by Env Psi content.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @return Data frame df updated with new column.
#' @export
compute_total_HIV_envPsi <- function(df){
  # extract ENV, Psi, und EnvPsi concentration from df (per Mio cells)
  env <- unique(df[grepl("Env", df$Target), "Mean Target/Mio cells"])
  psi <- unique(df[grepl("Psi", df$Target), "Mean Target/Mio cells"])
  envpsi <- unique(na.omit(df$`intact provirus/Mio cells Psi.Env, corrected for shearing`))
  df[["total HIV DNA/Mio cells (Env.Psi)"]] <- unlist(unname(env + psi - envpsi))
  return(df)
}

# create results tables (equivalent to analysis sheet in manual xlsx analysis)
#' Create analysis table for batch
#' 
#' Combine input data and compute analysis parameters for a group of wells (batch).
#' @param batch The wells which are analysed together.
#' @param conf_mat Confusion matrix (well vs. positives) for the batch.
#' @param num_target The total number of targets.
#' @param ch_dye Named list to match channel (names) with dyes (values).
#' @param multi_pos The possible multiple positives.
#' @param thresh The minimum number of accepted droplets in well to be included
#' in the analysis.
#' @param tar_mio_factor Factor to multiply Concentration with to obtain
#' Target/Mio cells
#' @param tab1 Data frame created in create_household_table().
#' @return List with two elements: a data frame and a sting indicating the nature
#' of the data sample (empty, H2O control, data)
#' @export
create_table <- function(batch, conf_mat, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1){
  
  # create table of respective wells
  tab <- merge(dtQC[dtQC$Well %in% batch, 1:8], conf_mat[conf_mat$Well %in% batch,], by="Well", all=T)
  tab <- Filter(function(x)!all(is.na(x)), tab)
  
  # get correct tar_mio_factor
  tar_mio <- tar_mio_factor[unique(tab$`Sample description 1`)]
  
  # change names (from 1_0_0_0 to Gag+ etc)
  markers <- match_channel_gene(tab, ch_dye, num_target)
  names(tab)[grepl("_", names(tab))] <- unname(sapply(grep("_", names(tab), value=T), transform_digits, markers))
  
  # do not compute values for H2O Wells
  if(unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")){
    return(list(tab, "h2o"))
  }
  
  # remove wells with less than x (7500) accepted droplets
  tab <- sufficient_droplets(tab, thresh)
    
  # check if any Well left after thresholding droplets
  if(nrow(tab) == 0){
    return(list("", "empty"))
  }
  
  # do not compute values for H2O Wells
  if(unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")){
    return(list(tab, "h2o"))
  }
  
  # add mean concentration and mean unsheared
  
  # check if sample description is in mean concentration household array
  if(unique(tab$`Sample description 1`) %in% tab1$`Sample description 1`){
    tab$`Mean concentration RPP30 + RPP30Shear (copies/µL)` <- unique(tab1[tab1$`Sample description 1` == unique(tab$`Sample description 1`), "Mean concentration RPP30 + RPP30Shear (copies/µL)"])[[1]]
      #mean_conc_household[unique(tab$`Sample description 1`)]
    tab$`Mean unsheared` <- unique(tab1[tab1$`Sample description 1` == unique(tab$`Sample description 1`), "Mean unsheared"])[[1]]
  } else{
    # if not, stop procedure
    stop("Sample description not exactly matched in names for Mean RPP30(Shear) concentration and Mean unsheared.
         Make sure to match Sample descriptions for RPP samples and HIV-target samples.")
  }

  # compute Target per million cells
  tab <- tab %>%
    mutate(`Target/Mio cells` = tar_mio * 10^6 * `Conc(copies/µL)`/(`Mean concentration RPP30 + RPP30Shear (copies/µL)`))
  
  # compute mean target per mio cells for each "group" (e.g. A09, B09, C09, D09)
  tab <- compute_groupwise_mean(tab, c("Target"), "Mean Target/Mio cells",
                                "Target/Mio cells")
  
  # compute values for multiple positives for specified combinations
  for(multip in multi_pos){
    tab <- get_multi_pos(tab, multip, tar_mio)
  }
  
  tab <- compute_total_HIV(tab, multi_pos)

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
create_tables <- function(grouped_data, in_csv, ch_dye, multi_pos, thresh, tar_mio_factor, tab1){

  # initialize empty lists to save output
  output_tables <- list()
  h2o_tables <- list()
  conf_mats <-  list()
  
  # loop over groups and run functions
  for (group in apply(grouped_data[is.na(grouped_data$RPP30), "shared_Wells"],
                      1,
                      function(x)(unname(unlist(strsplit(x, ", ")))))){
    # get relevant section of csv
    sub_in_csv <- in_csv[in_csv$Well %in% group,]
    # remove columns with all na
    sub_in_csv <- sub_in_csv[, !apply(is.na(sub_in_csv), 2, all)]
    # compute confusion matrix 
    # TODO: before conf_mat was computed once for all wells with the same targets
    # TODO: append appropriate conf matrices at the end
    conf_mat <- pivot_wider(sub_in_csv, id_cols = Well,
                            names_from = grep("Target", names(sub_in_csv), value = T),
                            values_from = Count,
                            values_fill = 0)

    # get number of targets
    num_target <- length(grep("Target", names(sub_in_csv), value = T))
    # check number of targets
    if(num_target < 4){
      warning("Less than 4 Targets detected. Skipping. If this is expexted 
              behaviour, please contact Mark :)")
      next
    }
    
    # compute output tables (equivalent to analysis sheet in manual analysis)

    # create output table
    out_table <- create_table(group, conf_mat, num_target, ch_dye, multi_pos, thresh, tar_mio_factor, tab1)
    
    # save to water table if water control
    if(out_table[[2]] == "h2o"){
      h2o_tables[[length(h2o_tables) + 1]] <- out_table[[1]]

    } else if(out_table[[2]] == "empty"){
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

# ======================== Beauty Output =======================================

#' Create analysis tables
#' 
#' Create an output sheet that collects the most important information.
#' @param table Data frame computed in create_tables().
#' @param multi_pos The possible multiple positives.
#' @return A named vector as blueprint for a xlsx sheet.
#' @export
get_output_sheet <- function(table, multi_pos, grouped_data){

  # create output of single positives
  sample <- unique(table$`Sample description 1`)
  output_sheet = c("Title" = paste(unique(table$Well), collapse = ", "),
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
                                             unique(table[["intact provirus/Mio cells Psi.Env, corrected for shearing"]])),
                   "Env-defective" = round(unique(table[grepl("Env", table$Target), "Mean Target/Mio cells"]) - 
                                             unique(table[["intact provirus/Mio cells Psi.Env, corrected for shearing"]]))
                   )
  
  # loop over multiple positives and add to output
  for (multip in multi_pos){
    name0 <- paste0("Intact concentration ", paste(multip, collapse = "."), " (copies/µl)")
    name1 <- paste0("intact provirus/Mio cells ", paste(multip, collapse = "."))
    name2 <- paste0(name1, ", corrected for shearing")

    x0 <- round(unique(na.omit(table[[name0]])), 2)
    x1 <- round(unique(na.omit(table[[name1]])))
    x2 <- round(unique(na.omit(table[[name2]])))

    output_sheet[[name0]] <- ifelse(length(x0) != 0,
                                    x0,
                                    NA)
    output_sheet[[name1]] <- ifelse(length(x1) != 0,
                                    x1,
                                    NA)
    output_sheet[[name2]] <- ifelse(length(x2) != 0,
                                    x2,
                                    NA)
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
write_output_file <- function(output_tables, conf_mats, tab1, output_file, h2o_table, multi_pos){
  l <- lapply(output_tables, get_output_sheet, multi_pos, grouped_data)
  output_sheet <- do.call(rbind.data.frame, l)
  
  
  output_list <- append(append(list(output_sheet), conf_mats), output_tables)
  output_list <- append(append(output_list, list(tab1)), h2o_table)

  openxlsx::write.xlsx(output_list, output_file)  
}



