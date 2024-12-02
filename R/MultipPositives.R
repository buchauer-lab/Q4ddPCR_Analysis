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
get_multi_pos <- function(df, genes, tar_mio_factor){
  # check that genes are specified
  if(is.null(genes)){
      stop("Need to specify genes to compute multiple positives.")
  }
  
  # create pattern to ask whether string contains all specified genes in arbitrary order
  pattern <- paste(paste0("(?=.*", genes, "\\+)"), collapse = "")

  # extract all column names of df that contain specified genes
  multi_pos_names <- grep(pattern, names(df), value = TRUE, perl = TRUE)

  # get counts positive for target
  tar_pos <- unlist(apply(df, 1, sum_target_positive_values))
  
  # get counts for multiple positives
  multi_pos <- apply(df[, multi_pos_names], 1, sum)
  
  # compute concentration
  conc_pos <- df$`Conc(copies/µL)` * as.vector(unlist(multi_pos)) / tar_pos
  
  # write concentration into df
  name <- paste0("Concentration ", paste(genes, collapse = "."),
                 " positive for target (copies/µl)")
  
  if(length(conc_pos) != 0){
    df[[name]] <- conc_pos
  } else{
    df[[name]] <- 0
  }
  
  # set NA to zero and not used Wells (Target not in multiple positive included) to NaN
  rows_w_target <- Reduce("|", lapply(genes, function(x)(grepl(x, df$Target))))
  df[is.na(df[, name]), name] = 0
  df[!(rows_w_target), name] = NaN
  
  # compute mean concentration by target and write into df
  df <- compute_groupwise_mean(df, c("Target"),
                                paste0("Mean ", name),
                                name)
  
  # compute mean by (target and) well and write into df
  name2 <- paste0("Intact concentration ", paste(genes, collapse = "."), " (copies/µl)")
  df <- compute_groupwise_mean(df, c(), name2,
                                paste0("Mean ", name))
  
  # compute intact provirus per million cells (and shearing corrected) and write into df
  name3 <- paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))
  df[[name3]] <- tar_mio_factor * 10^6 * df[[name2]] / tdfab$`Mean concentration RPP30 + RPP30Shear (copies/µL)`
  df[!(rows_w_target), c(name2, name3)] = NA
  df[[paste0(name3, ", corrected for shearing")]] <- df[[name3]] / df$`Mean unsheared`
  
  # return extended df
  return(df)
}

