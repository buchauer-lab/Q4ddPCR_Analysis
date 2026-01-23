#' Create RPP30 table
#'
#' Create the RPP30 table to compute Shearing factor and mean RPP30(Shear) concentration.
#' @param df data frame containing information on wells with RPP30 and RPP30Shear
#' @param mean_copies_factor Number to multiply Mean concentration RPP30 + Shear with to compute mean copies/well
#' @param mean_cells_per_reac_factor Factor to multiply Mean copies/cell with to obtain Mean cells per reaction (named list)
#' @import dplyr
#' @import tidyr
#' @return Updated data frame
#' @export
compute_shearing_factor <- function(df, mean_copies_factor, mean_cells_per_reac_factor) {
  # make some sanity checks
  # concentration corrected by dilution factor
  df <- df %>%
    mutate(
      `Concentration RPP30 (corrected by dilutionfactor) (copies/uL)` =
        `Dilution Factor` * `Conc(copies/uL)`
    )

  # compute mean concentrations (grouped by target and sample description)
  df <- df %>%
    group_by(group_id, Target) %>%
    # note: this is mean concentration for RPP30 and RPP30Shear individually
    # RPP30 in name could be removed
    mutate(
      `Mean concentration RPP30 (corrected by dilutionfactor) (copies/uL)` =
        mean(`Concentration RPP30 (corrected by dilutionfactor) (copies/uL)`, na.rm = TRUE)
    ) %>%
    mutate(
      `STD concentration RPP30 (corrected by dilutionfactor) (copies/uL)` =
        sd(`Concentration RPP30 (corrected by dilutionfactor) (copies/uL)`, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # compute mean concentration of PRR30 and RPP30Shear combined
    group_by(`Sample description 1`) %>%
    mutate(
      `Mean concentration RPP30 + RPP30Shear (copies/uL)` =
        mean(`Concentration RPP30 (corrected by dilutionfactor) (copies/uL)`, na.rm = TRUE)
    ) %>%
    ungroup()

  # compute unsheared value
  df$unsheared <- df$`Ch1+Ch2+` / (df$`Ch1+Ch2+` + (df$`Ch1+Ch2-` + df$`Ch1-Ch2+`) / 2)

  # compute mean unsheared
  df <- df %>%
    group_by(`Sample description 1`) %>%
    mutate(
      `Mean unsheared` =
        mean(`unsheared`, na.rm = TRUE)
    ) %>%
    mutate(
      `STD unsheared` =
        sd(`unsheared`, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Check if mean unsheared is 0 for some wells
  if(any(df$`Mean unsheared` == 0)){
    warning(paste0(
      "Mean unsheared is 0 for Wells: ", paste(unique(df$Well[df$`Mean unsheared` == 0]), collapse = ", "),
      ".
      Likely through all zeros for the Ch1+Ch2+ column in the input xlsx file. 
      The values will set to 1, i.e., no shearing correction will take place for
      the specified wells."
    ))
    df$`Mean unsheared`[df$`Mean unsheared` == 0] = 1
  }

  # compute shearing index
  df$`Shearing index` <- 1 - df$`Mean unsheared`

  # compute mean copies per well and mean cells per reaction
  df$`Mean copies/well` <- mean_copies_factor *
    df$`Mean concentration RPP30 + RPP30Shear (copies/uL)`
  df$`Mean cells per reaction` <- mean_cells_per_reac_factor[df$`Sample description 1`] *
    df$`Mean copies/well`

  return(df)
}

#' Define groups of well to analyse together
#'
#' Define groups based on the same sample description and targets found in this well
#' to analyse the group together.
#' @param dtQC output from read_files
#' @import dplyr
#' @import tidyr
#' @return Named list mapping well to group
#' @export
get_group_id <- function(dtQC){
  summary_df <- dtQC[, c("Well", "Sample description 1", "Target")] %>%
    group_by(Well) %>%
    summarise(
      description = first(`Sample description 1`),
      Target_set = paste(sort(unique(Target)), collapse = "|"),
      .groups = "drop"
    )
  id_df <- summary_df %>%
    group_by(description, Target_set) %>%
    mutate(group_id = cur_group_id())
  return(tibble::deframe(id_df[,c("Well", "group_id")]))
}

#' Add dilution factor
#'
#' Add the dilution factor to the table in define_groups. Gives warning if input objects do not match
#' @param df The data frame the dilution factor should be added on.
#' Must contain the column "Sample description 1".
#' @param dilution_factor A named vector that matches sample description (name)
#' with the used dilution factor (value).
#' @return Data frame updated with dilution factor
#' @export
add_dilution_factor <- function(df, dilution_factor) {
  # check if specified dilution factors match table
  if (any(!(df$`Sample description 1` %in% names(dilution_factor)))) {
    warning("Assuming 1 as dilution factor if not specified otherwise.")
    if (any(!(names(dilution_factor) %in% df$`Sample description 1`))) {
      warning("Not all specified dilution factor names appear in sample description.
                Please, check if this is expected.")
    }
    dilution_factor[unique(df$`Sample description 1`[!(df$`Sample description 1` %in% names(dilution_factor))])] <- 1
  }
  # add concentration
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
#' @export
sufficient_droplets <- function(df, thresh, droplet_col = "Accepted Droplets") {
  if (any(df[, droplet_col] < thresh)) {
    warning(paste0(
      "Removed well(s): ",
      paste(unique(df$`Well`[df[, droplet_col] < thresh]), collapse = ", "),
      ": Less than ", thresh, " droplets."
    ))
    df <- df[df[, droplet_col] >= thresh, ]
  }
  return(df)
}
