# branch main

# extract relevant information from csv and xlsx file

# load necessary packages
library(readxl)
library(tidyr)
library(dplyr)
source("MultipPositives.R")

# read files
read_files <- function(xlsx_file, csv_file, csv_skip, csv_nrows, remove_channel,
                       rm_zero_channels = FALSE){
  
  # check if files exist
  if(!file.exists(csv_file)){
    stop("Specified csv file does not exist. Please check if the path and name 
         are correct and that there are no spelling mistakes.")
  }
  
  if(!file.exists(xlsx_file)){
    stop("Specified xlsx file does not exist. Please check if the path and name 
         are correct and that there are no spelling mistakes.")
  }
  
  # read csv
  in_csv <- read.csv(csv_file, skip = csv_skip, row.names = NULL, 
                     nrows = csv_nrows)
  
  # remove shift in column names and resulting empty column at end
  colnames(in_csv) <- colnames(in_csv)[2:ncol(in_csv)]
  in_csv <- in_csv[,1:(ncol(in_csv) -1)]
  
  # remove wells as specified by user
  if(any(!(remove_channel %in% in_csv$Well))){
    stop("At least one of the Wells that should be removed does not occur in
         the Well column of the csv file. Please only remove Wells that do occur
         in the file.")
  }
  in_csv <- in_csv[!(in_csv$Well %in% remove_channel),]
  
  # read xlsx
  in_xlsx <- read_xlsx(xlsx_file, sheet = 1)
  
  # make concentration column numeric
  in_xlsx$`Conc(copies/µL)` <- as.numeric(in_xlsx$`Conc(copies/µL)`)
  in_xlsx$`Conc(copies/µL)`[is.na(in_xlsx$`Conc(copies/µL)`)] <- 0
  
  #remove wells as specified by user
  if(any(!(remove_channel %in% in_xlsx$Well))){
    stop("At least one of the Wells that should be removed does not occur in
         the Well column of the xlsx file. Please only remove Wells that do 
         occur in the file.")
  }
  in_xlsx <- in_xlsx[!(in_xlsx$Well %in% remove_channel),]
  
  # create sheet Data_table Quality Check (dtQC)
  cols_of_int <- c("Well", "Sample description 1", "DyeName(s)", "Target", 
                   "Conc(copies/µL)", "Accepted Droplets", "Positives", 
                   "Negatives", grep("Ch", names(in_xlsx), value = T))
    dtQC <- in_xlsx[,cols_of_int]
  
  # add threshold and number of total positives
  dtQC$Threshold <- dtQC$`Accepted Droplets`/3
  dtQC$`Total positives` <- apply(dtQC[,grep("\\+", names(dtQC), value=T)], 1, sum)
  
  # remove wells that have concentration zero for at least one well if specified
  if(rm_zero_channels){
    # make exception for H20 channels
    zero_ch <- unique(dtQC$Well[(dtQC$`Conc(copies/µL)` == 0) & 
                                  !(dtQC$`Sample description 1` %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser"))])
    if(length(zero_ch) > 0){
      warning(paste0("These wells will be removed, as they have at least one channel with concentration value 0: ",
                     paste(zero_ch, collapse = ", ")))
    }
    dtQC <- dtQC[!(dtQC$Well %in% zero_ch),]
    in_csv <- in_csv[!(in_csv$Well %in% zero_ch),]
  }
  
  # return read data
  return(list(in_csv, dtQC))
}

# create table with RPP30 and RPP30Shear controls
create_household_table <- function(dtQC, dilution_factor, custom_dilution_factor,
                                   thresh, mean_copies_factor, mean_cells_per_reac_factor){
  
  # check which wells have same staining
  table_targ_well<- table(dtQC[, c("Well", "Target")])
  df_tw <- as.data.frame.matrix(table_targ_well)
  print(head(df_tw))
  # get names of rows with Target RPP30 or RPP30Shear
  rpp30_names <- rownames(df_tw[which(df_tw$RPP30 == 1),])
  rpp30sh_names <- rownames(df_tw[which(df_tw$RPP30Shear == 1),])
  print(rpp30_names)
  print(rpp30sh_names)
  # check that wells for RPP30 and RPP30shear are the same
  if(any(rpp30_names != rpp30sh_names)){
    stop("Wells with targets RPP30 and RPP30Shear do not match. Please check
         again. If this is expected behaviour, please contact Mark :)")
  }
  
  # compute RPP30(Shear) table
  tab1 <- arrange(dtQC[dtQC$Well %in% rpp30_names,], Target)

  # remove Wells with less than x (7500) accepted droplets
  if(any(tab1$`Accepted Droplets` < thresh)){
    warning(paste0("Removed well(s): ", 
                 unique(tab1$`Accepted Droplets`[tab1$`Accepted Droplets` < thresh]),
                 ". Less than ", thresh, " droplets."))
    tab1 <- tab1[tab1$`Accepted Droplets` >= thresh,]
  }
  
  # add dilution factor
  if(custom_dilution_factor){
    # check if specified dilution factors match table
    if(any(!(tab1$`Sample description 1` %in% names(dilution_factor)))){
      warning("Assuming 1 as dilution factor if not specified otherwise.")
      if(any(!(names(dilution_factor) %in% tab1$`Sample description 1` ))){
        warning("Not all specified dilution factor names appear in sample description.")
      }
      
      dilution_factor[unique(tab1$`Sample description 1`[!(tab1$`Sample description 1` %in% names(dilution_factor))])] <- 1
    }
    
    # compute concentration
    tab1$`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)` <-
      dilution_factor[tab1$`Sample description 1`] * tab1$`Conc(copies/µL)` 
  
  } else {
    # use 1 as dilution factor for all concentrations
    warning("No dilution factors specified. Assuming 1 to be dilution factor.")
    tab1$`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)` <- tab1$`Conc(copies/µL)` 
  }
  
  # compute mean concentrations (grouped by target and sample description)
  print(tab1$`Sample description 1` %>% table)
  tab1 <- tab1 %>%
    group_by(Target, `Sample description 1`) %>%
    mutate(`Mean concentration RPP30 (corrected by dilutionfactor) (copies/µL)` =
             mean(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)`, na.rm = TRUE)) %>%
    mutate(`STD concentration RPP30 (corrected by dilutionfactor) (copies/µL)` =
             sd(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)`, na.rm = TRUE)) %>%
    ungroup() %>%
    # compute mean concentration of PRR30(Shear) (grouped by sample description)
    group_by(`Sample description 1`) %>%
    mutate(`Mean concentration RPP30 + RPP30Shear (copies/µL)` =
             mean(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)`, na.rm = TRUE)) %>%
    ungroup() 
  
  # compute unsheared value  
  tab1$unsheared <- tab1$`Ch1+Ch2+`/(tab1$`Ch1+Ch2+` + (tab1$`Ch1+Ch2-` + tab1$`Ch1-Ch2+`)/2)
  
  # compute mean unsheared
  tab1 <- tab1 %>%
    group_by(`Sample description 1`) %>%
    mutate(`Mean unsheared` =
             mean(`unsheared`, na.rm = TRUE)) %>%
    mutate(`STD unsheared` =
             sd(`unsheared`, na.rm = TRUE)) %>%
    ungroup() 
  
  # compute shearing index
  tab1$`Shearing index` <- 1-tab1$`Mean unsheared`
  
  # compute mean copies per well and mean cells per reaction
  tab1$`Mean copies/well` <- mean_copies_factor * tab1$`Mean concentration RPP30 + RPP30Shear (copies/µL)` # mean copies factor should be volume probably
  tab1$`Mean cells per reaction` <- mean_cells_per_reac_factor * tab1$`Mean copies/well`
  
  
  # extract mean concentration of RPP30 + RPP30Shear, mean unsheared, mean cells per reaction, and shearing index
  mean_conc_household <- setNames(unique(tab1$`Mean concentration RPP30 + RPP30Shear (copies/µL)`),
                                  unique(tab1$`Sample description 1`))
  mean_unsheared <- setNames(unique(tab1$`Mean unsheared`),
                             unique(tab1$`Sample description 1`))
  mean_cells_per_reac <- setNames(unique(tab1$`Mean cells per reaction`),
                                  unique(tab1$`Sample description 1`))
  Shear_ind <- setNames(unique(tab1$`Shearing index`),
                        unique(tab1$`Sample description 1`))
  
  return(list(tab1, mean_conc_household, mean_unsheared, mean_cells_per_reac, Shear_ind, rpp30_names, df_tw))
}


# ====================== Define functions to calculate table ==================

# function to get all possible combinations of multiple positives 
get_multipos <- function(markers){
  multi_pos <- list()
  
  for (i in 2:length(markers)) {
    # Get the combinations of length i
    comb <- combn(markers, i, simplify = FALSE)
    # Append the combinations to the list
    multi_pos <- c(multi_pos, comb)
  }
  
  return(multi_pos)
}

# helper function to convert the binary digits into their meaning
transform_digits <- function(combination, marker) {
  
  # Split combination into individual digits
  digits <- strsplit(combination, "_")[[1]]

  # Initialize empty vector to store object names
  positive_objects <- character(0)
  
  # Iterate over each digit and add corresponding object name if positive
  for (i in 1:length(digits)) {
    if (as.integer(digits[i]) == 1) {
      positive_objects <- c(positive_objects, marker[i])
    }
  }
  
  # Concatenate object names
  transformed_combination <- paste(positive_objects, collapse = "")
  if(transformed_combination == ""){
    return("All negative")
  }
  return(transformed_combination)
}

# helper function to extract colnames of target positive cols
extract_target_positive_colnames <- function(target, colnames){
  cols <- grep(target, colnames, value=T)
  return(grep(" ", cols, value=T, invert = T))
}

# helper function to get Channel - Gene match
match_channel_gene <- function(tab, ch_dye, num_ch){
  genes <- c()
  for (dye in ch_dye){
    gene <- tab[tab$"DyeName(s)" == dye, "Target"]
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

# helper function to compute total marker positive values
sum_target_positive_values <- function(tab){
  target <- substr(tab["Target"], (nchar(tab["Target"])-2), nchar(tab["Target"]))
  cols <- extract_target_positive_colnames(target, names(tab))
  return(sum(as.numeric(tab[cols]), na.rm=T))
}

# helper function to compute group specific mean (new_col) of a variable (old_col)
# e.g., mean of Target per Mio cells for A11 & B11 with ENV Target
compute_groupwise_mean<- function(tab, group_by, new_col, old_col){
  
  if(grepl("Mean", new_col)){
    tab <- tab %>%
      group_by_at(vars(!!!group_by)) %>%
      mutate(!!sym(new_col) :=
               mean(!!sym(old_col), na.rm = TRUE)) %>%
      mutate(!!sym(sub("Mean", "STD", new_col)) :=
               sd(!!sym(old_col), na.rm = TRUE)) %>%
      ungroup()
  } else{
    tab <- tab %>%
      group_by_at(vars(!!!group_by)) %>%
      mutate(!!sym(new_col) :=
               mean(!!sym(old_col), na.rm = TRUE)) %>%
      ungroup()
  }

  return(tab)
}

# Function to check if one combination is a subset of another
is_subset <- function(subset, set) {
  all(subset %in% set)
}

# Function to compute relevant subsets
get_subset <- function(multi_pos){
  # Initialize a list to store the subset relationships
  subset_list <- vector("list", length(multi_pos))
  names(subset_list) <- seq_along(multi_pos) 
  
  # Check which combinations are subsets of which other combinations
  for (i in 1:(length(multi_pos)-1)) {
    subset_list[[i]] <- list()
    for (j in 1:(length(multi_pos)-1)) {
      if (i != j) {
        if (is_subset(multi_pos[[i]], multi_pos[[j]])) {
          subset_list[[i]] <- c(subset_list[[i]], j)
        }
      }
    }
  }
  subset_list <- lapply(subset_list, unlist)
  return(subset_list)
}

# function to compute Total HIV DNA/Mio cells as proposed by Rachel
compute_total_HIV <- function(tab, multi_pos){
  # get length of multiple positives
  len_pos <- unlist(lapply(multi_pos, length))
  
  # find which multiple positives are a subset of another one
  subsets <- get_subset(multi_pos)

  # check that there are multiple positives that need to be corrected
  if(max(len_pos) < 2){
    warning("Maximum number of genes in multiple positives is 2. Total HIV DNA
            per Mio cells will not be computed")
    return(tab)
  }
  
  # get shear corrected all positive value
  quad_value_shear_corrected <- tab[[paste0("intact provirus/Mio cells ",
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
      original_value <- tab[[paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))]]
      # correct by all positive value
      new_value = original_value - quad_value_shear_corrected
      # correct by subsets if necessary
      if(!is.null(subsets[[as.character(comb)]])){
        for(subset in subsets[[as.character(comb)]]){
          # substract subset value
          subset_value <- tab[[paste0("Higher multipos corrected ", paste(multi_pos[[subset]], collapse = "."))]]
          #print(subset_value)
          new_value = new_value - subset_value
        }
      }
      # save new values to table
      tmp <- unique(na.omit(new_value))
      if(length(tmp) == 0){
        tmp = 0
      }
      tab[[name]] <- tmp
    }
    
  }
  # correct for single positives
  target <- substr(tab$Target, nchar(tab$Target)-2, nchar(tab$Target))

  # compute for each column
  subtract_multipos <- lapply(unique(target), function(tar){
    cols.of.int <-  grep(tar, grep("Higher multipos corrected", names(tab), value=T), value=T)
    res <- apply(tab[grepl(tar, tab$Target),cols.of.int], 1, sum, na.rm=TRUE)
    return(res)
  })

  tab[["Multipos corrected Mean Target/Mio cells"]] <- tab[["Mean Target/Mio cells"]] - quad_value_shear_corrected - unlist(subtract_multipos)
  
  tab[["Higher multipos corrected sum for all targets/Mio cells"]] <- 4 * mean(tab[["Multipos corrected Mean Target/Mio cells"]], na.rm=TRUE)
  
  # TODO
  tab[["total HIV DNA/Mio cells"]] <-  quad_value_shear_corrected + apply(tab[,grep("Higher multipos corrected", names(tab), value=TRUE)], 1, sum, na.rm=T)
  return(tab)
}

# create results tables (equivalent to analysis sheet in manual xlsx analysis)
create_table <- function(batch, conf_mat, num_target, ch_dye, multi_pos, thresh, tar_mio_factor){
  
  # create table of respective wells
  tab <- merge(dtQC[dtQC$Well %in% batch, 1:8], conf_mat[conf_mat$Well %in% batch,], by="Well", all=T)
  tab <- Filter(function(x)!all(is.na(x)), tab)
  
  # change names (from 1_0_0_0 to Gag+ etc)
  markers <- match_channel_gene(tab, ch_dye, num_target)
  names(tab)[grepl("_", names(tab))] <- unname(sapply(grep("_", names(tab), value=T), transform_digits, markers))
  
  # do not compute values for H2O Wells
  if(unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")){
    return(list(tab, "h2o"))
  }
  
  # remove wells with less than x (7500) accepted droplets
  if(any(tab$`Accepted Droplets` < thresh)){
    warning(paste0("Removed well: ", 
                   unique(tab$`Well`[tab$`Accepted Droplets` < thresh]),
                   ". Less than ", thresh, " droplets. "))
    
    tab <- tab[tab$`Accepted Droplets` >= thresh,]
    
    # check if any Well left
    if(nrow(tab) == 0){
      return(list("", "empty"))
    }
  }
  
  
  # do not compute values for H2O Wells
  if(unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")){
    return(list(tab, "h2o"))
  }
  
  # add mean concentration and mean unsheared
  
  # check if names of mean concentration and mean unsheared match
  if(!all(names(mean_conc_household) == names(mean_unsheared))){
    stop("Mean concentration of RPP30(Shear) and Mean unsheared are given for different Sample description.
         If this is expected behaviour, please contact Mark.")
  }
  
  # check if sample description is in mean concentration household array
  if(unique(tab$`Sample description 1`) %in% names(mean_conc_household)){
    mean_RPP <- mean_conc_household[unique(tab$`Sample description 1`)]
    mean_unsh <- mean_unsheared[unique(tab$`Sample description 1`)]
  } else{
    # if not, check if there is a partial match in names
    warning("Sample description not exactly matched in names for Mean RPP30(Shear) concentration and Mean unsheared.
            Trying flexibel appraoch next.")
    flex_name_bool <- unlist(lapply(names(mean_conc_household), grepl, unique(tab$`Sample description 1`)))
    print(names(mean_conc_household))
    print(unique(tab$`Sample description 1`))
    print(flex_name_bool)
    if(sum(flex_name_bool) != 1){
      # abort if there is no partial match and ask user to fix data
      stop("Flexibel approach failed. Make sure to match Sample descriptions for RPP samples and HIV-target samples.")
    }
    # print warning to show the partial match
    warning(paste0("Matched Mean conc RPP/Mean unsheared name ", names(mean_conc_household[flex_name_bool]),
                   " with sample description ", unique(tab$`Sample description 1`)))
    mean_RPP <- mean_conc_household[flex_name_bool]
    mean_unsh <- mean_unsheared[flex_name_bool]
  }
  
  # add mean concentration and mean unsheared to table
  tab$`Mean concentration RPP30 + RPP30Shear (copies/µL)` <- mean_RPP
  tab$`Mean unsheared` <- mean_unsh
  
  # compute Target per million cells
  tab <- tab %>%
    mutate(`Target/Mio cells` = tar_mio_factor * 10^6 * `Conc(copies/µL)`/(`Mean concentration RPP30 + RPP30Shear (copies/µL)`))
  
  # add temporary column with Well number (A09 -> 09)
  tab <- arrange(tab, by=Target, Well)
 
  
  # compute mean target per mio cells for each "group" (e.g. A09, B09, C09, D09)
  tab <- compute_groupwise_mean(tab, c("Target"), "Mean Target/Mio cells",
                                "Target/Mio cells")
  

  # compute values for multiple positives for specified combinations
  for(multip in multi_pos){
    tab <- get_multi_pos(tab, multip, tar_mio_factor)
  }
  
  tab <- compute_total_HIV(tab, multi_pos)

  return(list(tab, "tab"))
}



# ====================== Compute other stuff ===================================
define_groups <- function(df_tw, rpp30_names){
  # continue with other names
  df_tw <- df_tw[!(rownames(df_tw) %in% rpp30_names),]

  # add sample description
  table_samp_desc <- as.data.frame.matrix(table(dtQC[, c("Well" , "Sample description 1")]))
  comb_table <- cbind(df_tw, table_samp_desc[!(rownames(table_samp_desc) %in% rpp30_names),])
  
  # get wells that are used for same confusion matrix
  mat_sum <- apply(df_tw, 1, sum)
  mat_groups <- lapply(unique(mat_sum), function(x) names(which(mat_sum == x)))
  
  # get wells that should be analysed together
  comb_str <- apply(comb_table, 1, function(row) paste(row, collapse = "_"))
  groups <- lapply(unique(comb_str), function(x) names(which(comb_str == x)))
  
  # num_groups <- length(groups)
  return(list(mat_groups, groups))
}

# ====================== Execute ===============================================
create_tables <- function(mat_groups, in_csv, groups, ch_dye, multi_pos, thresh, tar_mio_factor){

  # initialize empty lists to save output
  output_tables <- list()
  h2o_tables <- list()
  conf_mats <-  list()
  
  # loop over groups and run functions
  for (mat_group in mat_groups){
    sub_in_csv <- in_csv[in_csv$Well %in% mat_group,]
    sub_in_csv <- sub_in_csv[, !apply(is.na(sub_in_csv), 2, all)]
    # compute confusion matrix (for 4 or 5 target group)
    conf_mat <- pivot_wider(sub_in_csv, id_cols = Well,
                            names_from = grep("Target", names(sub_in_csv), value = T),
                            values_from = Count)

    # set na to zero
    conf_mat[is.na(conf_mat)] <- 0
    
    # get number of targets
    num_target <- length(grep("Target", names(sub_in_csv), value = T))
    
    # check number of targets
    if(num_target < 4){
      warning("Less than 4 Targets detected. Skipping. If this is expexted 
              behaviour, please contact Mark :)")
      next
    }
    
    # compute output tables (equivalent to analysis sheet in manual analysis)
    for(group in groups){
      if(all(group %in% mat_group)){
        # create output table
        out_table <- create_table(group, conf_mat, num_target, ch_dye, multi_pos, thresh, tar_mio_factor)
        
        # save to water table if water control
        if(out_table[[2]] == "h2o"){
          h2o_tables[[length(h2o_tables) + 1]] <- out_table[[1]]

        } else if(out_table[[2]] == "empty"){
          # skip if empty (everything removed during QC)
            next
          } else {
          output_tables[[length(output_tables) + 1]] <- out_table[[1]]
        }
      } else if(all(!(group %in% mat_group))){
        next
      } else{
        stop("Not all samples in one confusion matrix. This is an error in the code! Contact Mark.")
      }
    }
    # save confusion matrix
    conf_mats[[length(conf_mats) + 1]] <- conf_mat
  }
  return(list(output_tables, conf_mats, h2o_tables))
}

# ======================== Beauty Output =======================================

# create output with the key results (results sheet in manual analysis)
get_output_sheet <- function(table, multi_pos){
  
  # check if mean cells per reaction name matches sample description or if partial match is necessary
  if(unique(table$`Sample description 1`) %in% names(mean_cells_per_reac)){
    desc <- names(mean_cells_per_reac[unique(table$`Sample description 1`)])
  } else{
    warning("Sample description not exactly matched in names for Mean RPP30(Shear) concentration and Mean unsheared.
            Trying flexibel appraoch next.")
    flex_name_bool <- unlist(lapply(names(mean_cells_per_reac), grepl, unique(table$`Sample description 1`)))
    if(sum(flex_name_bool) != 1){
      stop("Flexibel approach failed. Make sure to match Sample descriptions for RPP samples and HIV-target samples.")
    }
    warning(paste0("Matched Mean cells per reac name ", names(mean_cells_per_reac[flex_name_bool]),
                   " with sample description ", unique(table$`Sample description 1`)))
    desc <- names(mean_cells_per_reac[flex_name_bool])
  }
  # create output of single positives
  output_sheet = c("Title" = paste(unique(table$Well), collapse = ", "),
                   "Concentration probes" = unique(table$`Sample description 1`),
                   "Number of cells analysed" = unname(round(mean_cells_per_reac[desc])),
                   "Shearing Index" = unname(round(Shear_ind[desc], 2)),
                   "Total HIV-DNA (copies/Mio cells)" = round(unique(table[["total HIV DNA/Mio cells"]])),
                   "Gag/Mio cells" = round(unique(table[grepl("Gag", table$Target), "Mean Target/Mio cells"])),
                   "Pol/Mio cells" = round(unique(table[grepl("Pol", table$Target), "Mean Target/Mio cells"])),
                   "Psi/Mio cells" = round(unique(table[grepl("Psi", table$Target), "Mean Target/Mio cells"])),
                   "Env/Mio cells" = round(unique(table[grepl("Env", table$Target), "Mean Target/Mio cells"])),
                   "RU5/Mio cells" = ifelse(any(grepl("RU5", table$Target)),
                                            round(unique(table[grepl("RU5", table$Target), "Mean Target/Mio cells"])),
                                            NA)
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

# write information to xlsx file
write_output_file <- function(output_tables, conf_mats, tab1, output_file, h2o_table, multi_pos){
  l <- lapply(output_tables, get_output_sheet, multi_pos)
  output_sheet <- do.call(rbind.data.frame, l)
  
  
  output_list <- append(append(list(output_sheet), conf_mats), output_tables)
  output_list <- append(append(output_list, list(tab1)), h2o_table)
  
  openxlsx::write.xlsx(output_list, output_file)  
}



