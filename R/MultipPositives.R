#' Get count for a multiplet
#'
#' This function identifies the column name(s) of a multiplet based on a target 
#' list and returns the column values or the sum thereof
#' @param df data_table
#' @param genes list for which multiplets are computed for
#' @return (Sum of) multiplet count
#' @keywords internal
get_multiplet_count <- function(df, genes){
  # create pattern to ask whether string contains all specified genes in arbitrary order
  pattern <- paste(paste0("(?=.*", genes, "\\+)"), collapse = "")
  # extract all column names of df that contain specified genes
  multi_pos_names <- grep(pattern, names(df), value = TRUE, perl = TRUE)
  
  # get counts for multiple positives
  if(length(multi_pos_names) == 1){
    multi_pos <- df[, multi_pos_names]
  } else if (length(multi_pos_names) == 0){
    warning(paste("No multiple positives were found in data for:"), genes,
            " in well(s): ", paste(unique(df$Well), collapse=", "))
    multi_pos <- 0
  } else{
    multi_pos <- apply(df[, multi_pos_names], 1, sum)
  }
  return(multi_pos)
}

# Function to return concentration of multiple positive samples (e.g. Env+Psi+ or quadrupel positive)
#' Get multiple positives concentrations
#'
#' Compute the (intact) concentration and its mean across grouped Well(s) for
#' multiple positives (e.g., Env+Psi+, Env+Gag+Pol+Psi+).
#' @param df Data frame with the information on the measurement.
#' @param genes The genes for which the sample should be positive
#' (e.g., c("Env", "Psi") for Env+Psi+
#' @param tar_mio_factor Factor to multiply Concentration with to obtain
#'  Target/Mio cells
#' @return Updated data frame (df)
#' @export
get_multi_pos <- function(df, genes, tar_mio_factor) {
  # check that genes are specified
  if (is.null(genes)) {
    stop("Need to specify genes to compute multiple positives.")
  }
  # get mutli_pos
  multi_pos <- get_multiplet_count(df, genes)
  # compute concentration
  # df$Positives: get counts positive for target (row-specific)
  conc_pos <- df$`Conc(copies/uL)` * as.vector(unlist(multi_pos)) / df$Positives
  
  # write concentration into df
  name <- paste0(
    "Concentration ", paste(genes, collapse = "."),
    " positive for target (copies/ul)"
  )

  if (length(conc_pos) != 0) {
    df[[name]] <- conc_pos
  } else {
    df[[name]] <- 0
  }
  # set NA to zero and not used Wells (Target not in multiple positive included) to NaN
  rows_w_target <- Reduce("|", lapply(genes, function(x) (grepl(x, df$Target))))
  df[is.na(df[, name]), name] <- 0
  df[!(rows_w_target), name] <- NA

  # compute mean concentration by target and write into df
  df <- df %>%
    group_by(group_id, Target) %>%
    mutate(!!sym(paste0('Mean ', name)) := mean(!!sym(name), na.rm = TRUE)) %>%
    mutate(!!sym(paste0('SD ', name)) := sd(!!sym(name), na.rm = TRUE)) %>%
    ungroup()
    

  # compute mean per group and write into df
  name2 <- paste0("Intact concentration ", paste(genes, collapse = "."), " (copies/ul)")
  df <- df %>%
    group_by(group_id) %>%
    mutate(!!sym(name2) := mean(!!sym(name), na.rm = TRUE)) %>%
    mutate(!!sym(paste0('SD ', name2)) := sd(!!sym(name), na.rm = TRUE)) %>%
    ungroup()

  # compute intact provirus per million cells (and shearing corrected) and write into df
  name3 <- paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))
  df <- df %>%
    mutate(!!sym(name3) := tar_mio_factor[`Sample description 1`] * 10^6 * !!sym(name2) / (`Mean concentration RPP30 + RPP30Shear (copies/uL)`)) %>%
    mutate(!!sym(paste0("SD ", name3)) := abs(tar_mio_factor[`Sample description 1`] * 10^6 / `Mean concentration RPP30 + RPP30Shear (copies/uL)`) * !!sym(paste0("SD ", name2)))
  
  df[!(rows_w_target), c(name2, name3)] <- NA
  
  df[[paste0(name3, ", corrected for shearing")]] <- df[[name3]] / df$`Mean unsheared`
  df[[paste0("SD ", name3, ", corrected for shearing")]] <- df[[paste0("SD ", name3)]] / df$`Mean unsheared`

  # return extended df
  return(df)
}

