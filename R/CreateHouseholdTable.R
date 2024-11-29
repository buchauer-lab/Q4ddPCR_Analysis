#' Create QC table
#' 
#' Create the dtQC table, which contains the relevant content from the xlsx file.
#' @param xlsx_df The data frame created when reading the xlsx file.
#' @param cols_of_int Columns that should be contained in dtQC. If not specified,
#' a standard set of columns is kept. (optional)
#' @param acc_drop_factor Numeric factor by which the threshold is computed from
#' the accepted droplets. (default=3)
#' @return The dtQC dataframe with the relevant columns.
#' @export
create_dtQC <- function(xlsx_df, cols_of_int = NULL, acc_drop_factor = 3){

}


# create table with RPP30 and RPP30Shear controls
create_household_table <- function(dtQC, grouped_data,
                                   thresh, mean_copies_factor, mean_cells_per_reac_factor){
  
  # compute RPP30(Shear) table
  tab1 <- arrange(dtQC[dtQC$Target %in% c("RPP30", "RPP30Shear"),], Target)
  
  # remove Wells with less than x (7500) accepted droplets
  tab1 <- sufficient_droplets(tab1, thresh)
  
  #  add dilution factor to tab1
  tab1 <- tab1 %>%
    left_join(grouped_data[!is.na(grouped_data$RPP30),] %>%
                select(`Sample description 1`, `Dilution Factor`),
              by = "Sample description 1")
  
  # compute concentration corrected by dilution factor
  tab1 <- tab1 %>%
    mutate(`Concentration RPP30 (corrected by dilutionfactor) (copies/µL)` =
             `Dilution Factor` * `Conc(copies/µL)`)
  
  # compute mean concentrations (grouped by target and sample description)
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
  tab1$`Mean copies/well` <- mean_copies_factor *
    tab1$`Mean concentration RPP30 + RPP30Shear (copies/µL)` 
  tab1$`Mean cells per reaction` <- mean_cells_per_reac_factor *
    tab1$`Mean copies/well`
  
  # add relevant information to grouped data
  grouped_data <- grouped_data %>%
    left_join(tab1 %>% select(`Sample description 1`,
                              `Mean concentration RPP30 + RPP30Shear (copies/µL)`,
                              `Mean unsheared`,
                              `Mean cells per reaction`,
                              `Shearing index`) %>%
                distinct(),
              by = "Sample description 1")
  
  
  return(list(tab1, grouped_data))
}


define_groups <- function(dtQC, dilution_factor){
  # make one column per Well
  transformed_data <- dtQC[, c("Well" , "Sample description 1", "Target")] %>%
    pivot_wider(names_from = Target, values_from = Target)
  # find wells with same sample description and targets
  grouped_data <- transformed_data %>%
    group_by(across(names(transformed_data)[2:ncol(transformed_data)])) %>%
    summarise(shared_Wells = paste(Well, collapse = ", "), .groups = "drop")
  
  # check if rpp30 and rpp30shear appear always together
  if (any(xor(is.na(grouped_data$RPP30), is.na(grouped_data$RPP30Shear)))){
    stop("Wells with targets RPP30 and RPP30Shear do not match. Please check
         again. If this is expected behaviour, please contact Mark :)")
  }
  # add the dilution factor if specified
  if(is_named_vector(dilution_factor)){
    grouped_data <- add_dilution_factor(grouped_data, dilution_factor)
  } else {
    warning("No dilution factor has been specified or not as named vector.
            Using dilution factor 1 for all samples.")
    grouped_data$`Dilution Factor` <- 1
  }
  return(grouped_data)
}

# check if an object is a named vector
is_named_vector <- function(obj) {
  # Check if the object is an atomic vector
  is_atomic_vector <- is.atomic(obj) && !is.null(obj)
  
  # Check if the object has non-empty names
  has_names <- !is.null(names(obj)) && all(nzchar(names(obj)))
  
  # Return TRUE only if both conditions are satisfied
  is_atomic_vector && has_names
}

# function to add the dilution factor to the table
add_dilution_factor <- function(tab, dilution_factor){
    # check if specified dilution factors match table
    if(any(!(tab$`Sample description 1` %in% names(dilution_factor)))){
      warning("Assuming 1 as dilution factor if not specified otherwise.")
      if(any(!(names(dilution_factor) %in% tab1$`Sample description 1` ))){
        warning("Not all specified dilution factor names appear in sample description.
                Please, check if this is expected.")
      }
      dilution_factor[unique(tab$`Sample description 1`[!(tab$`Sample description 1` %in% names(dilution_factor))])] <- 1
    }
    # compute concentration
    tab$`Dilution Factor` <- dilution_factor[tab$`Sample description 1`]

  return(tab)
}

# function to check if the number of droplets is above the minimum threshold
sufficient_droplets <- function(tab, thresh, droplet_col="Accepted Droplets"){
  if(any(tab[, droplet_col] < thresh)){
    warning(paste0("Removed well(s): ", 
                   paste(unique(tab$`Well`[tab[, droplet_col] < thresh]), collapse = ", "),
                   ": Less than ", thresh, " droplets."))
    tab <- tab[tab[, droplet_col] >= thresh,]
  }
  return(tab)
}