# branch dev
# extract relevant information from csv and xlsx file

# load necessary packages
library(readxl)
library(tidyr)
library(dplyr)

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
  # TODO: find more elegant solution, problem arises due to row.names = NULL
  colnames(in_csv) <- colnames(in_csv)[2:ncol(in_csv)]
  in_csv <- in_csv[,1:(ncol(in_csv) -1)]
  
  # remove wells with zero channels
  if(any(!(remove_channel %in% in_csv$Well))){
    stop("At least one of the Wells that should be removed does not occur in
         the Well column of the csv file. Please only remove Wells that do occur
         in the file.")
  }
  in_csv <- in_csv[!(in_csv$Well %in% remove_channel),]
  
  # read xlsx
  in_xlsx <- read_xlsx(xlsx_file, sheet = 1)
  
  # make Concentration col numeric
  in_xlsx$`Conc(copies/µL)` <- as.numeric(in_xlsx$`Conc(copies/µL)`)
  in_xlsx$`Conc(copies/µL)`[is.na(in_xlsx$`Conc(copies/µL)`)] <- 0
  
  #remove wells with zero channel
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
  
  # remove wells that have concentration zero for at least one well
  if(rm_zero_channels){
    zero_ch <- unique(dtQC$Well[(dtQC$`Conc(copies/µL)` == 0) & 
                                  !(dtQC$`Sample description 1` %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser"))])
    if(length(zero_ch) > 0){
      warning(paste0("These wells will be removed, as they have at least one channel with concentration value 0: ",
                     paste(zero_ch, collapse = ", ")))
    }
    dtQC <- dtQC[!(dtQC$Well %in% zero_ch),]
    in_csv <- in_csv[!(in_csv$Well %in% zero_ch),]
  }
  
  return(list(in_csv, dtQC))
}

# create table with RPP30 and RPP30Shear controls
create_household_table <- function(dtQC, dilution_factor, custom_dilution_factor){
  
  # check which wells have same staining
  table_targ_well<- table(dtQC[, c("Well", "Target")])
  
  df_tw <- as.data.frame.matrix(table_targ_well)
  
  #  Compute concentration of unsheared DNA 
  rpp30_names <- rownames(df_tw[which(df_tw$RPP30 == 1),])
  rpp30sh_names <- rownames(df_tw[which(df_tw$RPP30Shear == 1),])
  
  # check that wells for rpp30 and rpp30shear are the same
  if(any(rpp30_names != rpp30sh_names)){
    stop("Wells with targets RPP30 and RPP30Shear do not match. Please check
         again. If this is expected behaviour, please contact Mark :)")
  }
  
  # compute table 1 (RPP30 & RPP30Shear values)
  tab1 <- arrange(dtQC[dtQC$Well %in% rpp30_names,], Target)
  
  # remove Wells with less than 7500 accepted droplets
  if(any(tab1$`Accepted Droplets` < 7500)){
    warning(paste0("Removed well(s): ", 
                 unique(tab1$`Accepted Droplets`[tab1$`Accepted Droplets` < 7500]),
                 ". Less than 7500 droplets."))
    tab1 <- tab1[tab1$`Accepted Droplets` >= 7500,]
  }
  
  # add dilution factor
  if(custom_dilution_factor){
    # check if specified dilution factors match table
    if(any(!(tab1$`Sample description 1` %in% names(dilution_factor)))){
      warning("Assuming 1 as dilution factor if not specified otherwise.")
      d <- rep(1, length(dilution_factor))
      d[names(dilution_factor)] <- dilution_factor
      dilution_factor <- d
    }
    tab1$`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)` <-
      dilution_factor[tab1$`Sample description 1`] * tab1$`Conc(copies/µL)` 
  } else {
    warning("No dilution factors specified. Assuming 1 to be dilution factor.")
    tab1$`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)` <- tab1$`Conc(copies/µL)` 
  }
  
  # compute mean concentrations
  tab1 <- tab1 %>%
    group_by(Target, `Sample description 1`) %>%
    mutate(`Mean concentration RPP30 (corrected by dilutionfactor) (copies/µL)` =
             mean(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)`, na.rm = TRUE)) %>%
    mutate(`STD concentration RPP30 (corrected by dilutionfactor) (copies/µL)` =
             sd(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)`, na.rm = TRUE)) %>%
    ungroup() %>%
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
  tab1$`Mean copies/well` <- 20*tab1$`Mean concentration RPP30 + RPP30Shear (copies/µL)`
  tab1$`Mean cells per reaction` <- 2*tab1$`Mean copies/well`
  
  
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
  return(grep(target, colnames, value=T))
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
  if(new_col == "Mean Target/Mio cells"){

  }
  tab <- tab %>%
    group_by_at(vars(!!!group_by)) %>%
    mutate(!!sym(new_col) :=
             mean(!!sym(old_col), na.rm = TRUE)) %>%
    ungroup()
  return(tab)
}

# create results tables (equivalent to analysis sheet in manual xlsx analysis)
create_table <- function(batch, conf_mat, num_target, ch_dye){
  
  # create table of respective wells
  tab <- merge(dtQC[dtQC$Well %in% batch, 1:8], conf_mat[conf_mat$Well %in% batch,], by="Well", all=T)
  tab <- Filter(function(x)!all(is.na(x)), tab)

  # remove wells with less than 7500 accepted droplets
  if(any(tab$`Accepted Droplets` < 7500)){
    warning(paste0("Removed well: ", 
                   unique(tab$`Well`[tab$`Accepted Droplets` < 7500]),
                   ". Less than 7500 droplets. "))
    
    tab <- tab[tab$`Accepted Droplets` >= 7500,]
    
    # check if any Well left
    if(nrow(tab) == 0){
      return(list("", "empty"))
    }
  }
  
  # change names (from 1_0_0_0 to Gag+ etc)
  markers <- match_channel_gene(tab, ch_dye, num_target)
  names(tab)[grepl("_", names(tab))] <- unname(sapply(grep("_", names(tab), value=T), transform_digits, markers))
  
  # do not compute values for H2O Wells
  if(unique(tab$`Sample description 1`) %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser")){
    return(list(tab, "h2o"))
  }
  
  # add mean concentration and mean unsheared
  tab$`Mean concentration RPP30 + RPP30Shear (copies/µL)` <- mean_conc_household[unique(tab$`Sample description 1`)]
  tab$`Mean unsheared` <- mean_unsheared[unique(tab$`Sample description 1`)]
  
  # compute Target per million cells
  tab <- tab %>%
    mutate(`Target/Mio cells` = 2 * 10^6 * `Conc(copies/µL)`/(`Mean concentration RPP30 + RPP30Shear (copies/µL)`))
  
  # add temporary column with Well number (A09 -> 09)
  tab$tmp <- substr(tab$Well, 2, 3)
  tab <- arrange(tab, by=Target, tmp)
  
  # compute mean target per mio cells for each "group" (e.g. A09, B09, C09, D09)
  tab <- compute_groupwise_mean(tab, c("Target"), "Mean Target/Mio cells",
                                "Target/Mio cells")
  
  # compute all positive concentration (quadruple or quintuple)
  all_pos_name <- grep(paste(c(rep(".*\\+", num_target), ".*"), collapse = ""), names(tab), value = TRUE)
  tar_pos <- unlist(apply(tab, 1, sum_target_positive_values))
  all_pos <- tab[, all_pos_name]

  conc_all_pos <- tab$`Conc(copies/µL)` * as.vector(unlist(all_pos))/tar_pos
  if(length(conc_all_pos) != 0){
    tab$`Concentration all positive for target (copies/µL)` = conc_all_pos 
  } else {
    tab$`Concentration all positive for target (copies/µL)` = 0
  }
  
  
  # compute mean of all positive concentration
  tab <- compute_groupwise_mean(tab, c("Target"),
                                "Mean concentration all positive all targets (copies/µL)",
                                "Concentration all positive for target (copies/µL)")
  
  tab <- compute_groupwise_mean(tab, c(), "Intact concentration (copies/µl)",
                                "Mean concentration all positive all targets (copies/µL)")
  
  tab <- tab %>%
    mutate(`intact provirus/Mio cells` = 2 * 10^6 * `Intact concentration (copies/µl)` / `Mean concentration RPP30 + RPP30Shear (copies/µL)`)
  
  tab <- tab %>%
    mutate(`intact provirus/Mio cells, corrected for shearing` = `intact provirus/Mio cells` / `Mean unsheared`)
  
  # compute quadruple positives for 5 Target assays
  if(num_target == 5){
    quad_pos_names <- grep("(?=.*Gag\\+)(?=.*Pol\\+)(?=.*Psi\\+)(?=.*Env\\+)", names(tab), value = TRUE, perl = TRUE)
    tar_pos <- unlist(apply(tab, 1, sum_target_positive_values))
    quad_pos <- apply(tab[, quad_pos_names], 1, sum)

    conc_quad_pos <- tab$`Conc(copies/µL)` * as.vector(unlist(quad_pos))/tar_pos
    if(length(conc_all_pos) != 0){
      tab$`Concentration quad positive for target (copies/µL)` = conc_quad_pos 
    } else {
      tab$`Concentration quad positive for target (copies/µL)` = 0
    }
    
    # do not compute for RU5 wells
    tab[tab$Target == "RU5", "Concentration quad positive for target (copies/µL)"] = NA
    
    tab <- compute_groupwise_mean(tab, c("Target"),
                                  "Mean concentration quad positive all targets (copies/µL)",
                                  "Concentration quad positive for target (copies/µL)")
    
    tab <- compute_groupwise_mean(tab, c("tmp"), "Intact concentration quad (copies/µl)",
                                  "Mean concentration quad positive all targets (copies/µL)")
    
    tab <- tab %>%
      mutate(`intact provirus/Mio cells quad` = 2 * 10^6 * `Intact concentration quad (copies/µl)` / `Mean concentration RPP30 + RPP30Shear (copies/µL)`)
    
    tab <- tab %>%
      mutate(`intact provirus/Mio cells quad, corrected for shearing` = `intact provirus/Mio cells quad` / `Mean unsheared`)
    
  }
  
  # Psi and Env concentration
  Psip <- tab[,grep("Psi", names(tab), value=T)] %>% apply(1,sum, na.rm=T)
  Envp <- tab[,grep("Env", names(tab), value=T)] %>% apply(1,sum, na.rm=T)
  EnvpPsip <- tab[,grep("Env.*Psi|Psi.*Env", names(tab), value=T)] %>% apply(1,sum, na.rm=T)
  
  tab$`Concentration Env+Psi positive for target (copies/µL)` <- 
    ifelse(grepl("Psi", tab$Target),
           tab$`Conc(copies/µL)` * EnvpPsip / Psip, NA)
  
  tab$`Concentration Env+Psi positive for target (copies/µL)`[grepl("Env", tab$Target)] <- 
    tab$`Conc(copies/µL)`[grepl("Env", tab$Target)] * EnvpPsip[grepl("Env", tab$Target)] / Envp[grepl("Env", tab$Target)]
  
  tab$`Env+Psi+ concentration (copies/µL)` <- mean(tab$`Concentration Env+Psi positive for target (copies/µL)`, na.rm=T)
  
  tab$`Env+Psi+ concentration (copies/µL)`[is.na(tab$`Concentration Env+Psi positive for target (copies/µL)`)] <- NaN
  
  tab <- tab %>%
    mutate(`Env+Psi+ intact provirus/Mio cells` = 2 * 10^6 * `Env+Psi+ concentration (copies/µL)` / `Mean concentration RPP30 + RPP30Shear (copies/µL)`)
  
  tab <- tab %>%
    mutate(`Env+Psi+ intact provirus/Mio cells, Shear corrected` =  `Env+Psi+ intact provirus/Mio cells` / `Mean unsheared`)
  
  
  
  # remove tmp column
  tab <- tab[, -which(names(tab) == "tmp")]
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
create_tables <- function(mat_groups, in_csv, groups, ch_dye){
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
    conf_mat[is.na(conf_mat)] <- 0
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
        out_table <- create_table(group, conf_mat, num_target, ch_dye)
        if(out_table[[2]] == "h2o"){
          h2o_tables[[length(h2o_tables) + 1]] <- out_table[[1]]
        } else if(out_table[[2]] == "empty"){
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
    conf_mats[[length(conf_mats) + 1]] <- conf_mat
  }
  return(list(output_tables, conf_mats, h2o_tables))
}

# ======================== Beauty Output =======================================

# create output with the key results (results sheet in manual analysis)
get_output_sheet <- function(table){
  #desc <- sub("-.*", "", unique(table$`Sample description 1`)) 
  desc <- unique(table$`Sample description 1`)
  output_sheet = c("Title" = paste(unique(table$Well), collapse = ", "),
                   "Concentration probes" = unique(table$`Sample description 1`),
                   "Number of cells analysed" = unname(round(mean_cells_per_reac[desc])),
                   "Shearing Index" = unname(round(Shear_ind[desc], 2)),
                   "Total HIV-DNA (copies/Mio cells)" = ifelse(any("RU5" %in% unique(table$Target)),
                                                               round(unique(table[table$Target == "RU5", "Mean Target/Mio cells" ])),
                                                               NA),
                   "Concentration quadrupel positive (copies/µl)" = ifelse(any(grepl("quad", names(table))),
                                                                           round(unique(table$`Intact concentration quad (copies/µl)`), 2),
                                                                           round(unique(table$`Intact concentration (copies/µl)`), 2)),
                   "Concentration quintuple positive (copies/µl)" = ifelse(any(grepl("quad", names(table))),
                                                                           round(unique(table$`Intact concentration (copies/µl)`), 2),
                                                                           NA),
                   "Gag/Mio cells" = round(unique(table[grepl("Gag", table$Target), "Mean Target/Mio cells"])),
                   "Pol/Mio cells" = round(unique(table[grepl("Pol", table$Target), "Mean Target/Mio cells"])),
                   "Psi/Mio cells" = round(unique(table[grepl("Psi", table$Target), "Mean Target/Mio cells"])),
                   "Env/Mio cells" = round(unique(table[grepl("Env", table$Target), "Mean Target/Mio cells"])),
                   "RU5/Mio cells" = ifelse(any(grepl("RU5", table$Target)),
                                            round(unique(table[grepl("RU5", table$Target), "Mean Target/Mio cells"])),
                                            NA),
                   "Env+Psi+/Mio cells uncorr" = max(0, round(unique(na.omit(table$`Env+Psi+ intact provirus/Mio cells`)))), # will return 0 if value is NULL
                   "4c intact provirus/Mio cells uncorr" = ifelse(any(grepl("quad", names(table))),
                                                                  round(unique(table$`intact provirus/Mio cells quad`)),
                                                                  round(unique(table$`intact provirus/Mio cells`))),
                   "5c intact provirus/Mio cells uncorr" = ifelse(any(grepl("quad", names(table))),
                                                                  round(unique(table$`intact provirus/Mio cells`)),
                                                                  NA),
                   "Env+Psi+/Mio cells Shear corr" = max(0, round(unique(na.omit(table$`Env+Psi+ intact provirus/Mio cells, Shear corrected`)))), # will return 0 if value is NULL
                   "4c intact provirus/Mio cells Shear corr" = ifelse(any(grepl("quad", names(table))),
                                                                      round(unique(table$`intact provirus/Mio cells quad, corrected for shearing`)),
                                                                      round(unique(table$`intact provirus/Mio cells, corrected for shearing`))),
                   "5c intact provirus/Mio cells Shear corr" = ifelse(any(grepl("quad", names(table))),
                                                                      round(unique(table$`intact provirus/Mio cells, corrected for shearing`)),
                                                                      NA)
                   )
  return(output_sheet)
}

# =================== Save to xlsx file =======================================

# write information to xlsx file
write_output_file <- function(output_tables, conf_mats, tab1, output_file, h2o_table){
  l <- lapply(output_tables, get_output_sheet)
  output_sheet <- do.call(rbind.data.frame, l)
  
  
  output_list <- append(append(list(output_sheet), conf_mats), output_tables)
  output_list <- append(append(output_list, list(tab1)), h2o_table)
  
  openxlsx::write.xlsx(output_list, output_file)  
}



