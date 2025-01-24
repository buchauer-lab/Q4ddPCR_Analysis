#' Create RPP30 table
#' 
#' Create the RPP30 table to compute Shearing factor and mean RPP30(Shear) concentration.
#' @param dtQC The data frame created in the Readers/read_files function.
#' @param grouped_data Data frame created in define_groups function.
#' @param thresh Minimum number of accepted droplets to include the well in the calculations.
#' Other wise it will be removed during the process.
#' @param mean_copies_factor Number to multiply Mean concentration RPP30 + Shear with to compute mean copies/well # I guess should be named volume
#' @param mean_cells_per_reac_factor Factor to multiply Mean copies/cell with to obtain Mean cells per reaction
#' @import dplyr
#' @return List with two elements: data frame with information about RPP30(Shear) and an updated grouped_data.
#' @export
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

#' Define groups
#' 
#' Define the groups of Well that are analysed together. This is based upon the
#' columns "Sample description 1" and "Target". If they are the same for two
#' wells, they will be analysed together.
#' @param dtQC The data frame created in the Readers/read_files function.
#' @param dilution_factor The dilution factor will be applied to the concentration
#' to correct for dilution during the experiment.
#' @import dplyr
#' @return grouped_data data frame, with basic information to each found group.
#' @export
define_groups <- function(dtQC, dilution_factor){
  # make one column per Well
  transformed_data <- dtQC[, c("Well" , "Sample description 1", "Target")] %>%
    tidyr::pivot_wider(names_from = Target, values_from = Target)
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

#' Check if is named vector
#' 
#' Check if an object is a named vector
#' @param obj The object to check
#' @return Boolean
is_named_vector <- function(obj) {
  # Check if the object is an atomic vector
  is_atomic_vector <- is.atomic(obj) && !is.null(obj)
  
  # Check if the object has non-empty names
  has_names <- !is.null(names(obj)) && all(nzchar(names(obj)))
  
  # Return TRUE only if both conditions are satisfied
  is_atomic_vector && has_names
}

#' Add dilution factor
#' 
#' Add the dilution factor to the table in define_groups. Gives warning if input objects do not match
#' @param df The data frame the dilution factor should be added on.
#' Must contain the column "Sample description 1".
#' @param dilution_factor A named vector that matches sample description (name)
#' with the used dilution factor (value).
#' @return Data frame updated with dilution factor
add_dilution_factor <- function(df, dilution_factor){
    # check if specified dilution factors match table
    if(any(!(df$`Sample description 1` %in% names(dilution_factor)))){
      warning("Assuming 1 as dilution factor if not specified otherwise.")
      if(any(!(names(dilution_factor) %in% df$`Sample description 1` ))){
        warning("Not all specified dilution factor names appear in sample description.
                Please, check if this is expected.")
      }
      dilution_factor[unique(df$`Sample description 1`[!(df$`Sample description 1` %in% names(dilution_factor))])] <- 1
    }
    # compute concentration
    df$`Dilution Factor` <- dilution_factor[df$`Sample description 1`]

  return(df)
}

#' Check for sufficient droplets
#' 
#' Check if the number of accepted droplets is above the minimum threshold.
#' @param df The data frame that should be checked. Must contain the columns
#'  "Well" and the one specified in droplet_col.
#' @param thresh The threshold (minimum of accepted droplets).
#' @param droplet_col The column containing the information about the accepted
#' droplets (default: "Accepted Droplets")
#' @return Data frame with Wells above the given threshold only.
sufficient_droplets <- function(df, thresh, droplet_col="Accepted Droplets"){
  if(any(df[, droplet_col] < thresh)){
    warning(paste0("Removed well(s): ", 
                   paste(unique(df$`Well`[df[, droplet_col] < thresh]), collapse = ", "),
                   ": Less than ", thresh, " droplets."))
    df <- df[df[, droplet_col] >= thresh,]
  }
  return(df)
}