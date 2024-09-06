# Function to return concentration of multiple positive samples (e.g. Env+Psi+ or quadrupel positive)
get_multi_pos <- function(tab, genes, tar_mio_factor){
  # check that genes are speicfied
  if(is.null(genes)){
      stop("Need to specify genes to compute multiple positives.")
  }
  
  # create pattern to ask whether string contains all specified genes in arbitrary order
  pattern <- paste(paste0("(?=.*", genes, "\\+)"), collapse = "")

  # extract all column names of table that contain specified genes
  multi_pos_names <- grep(pattern, names(tab), value = TRUE, perl = TRUE)

  # Check that only count columns are in, not previously computed multipos concentrations
  # e.g., do not use "Env+Psi+Gag+Pol+ Intact Provirus" when computing "Env+Psi+ concentration"
  multi_pos_names <- grep(" ", multi_pos_names, value = TRUE, invert = TRUE)
  
  # get counts positive for target
  tar_pos <- unlist(apply(tab, 1, sum_target_positive_values))
  
  # get counts for multiple positives
  multi_pos <- apply(tab[, multi_pos_names], 1, sum)
  
  # compute concentration
  conc_pos <- tab$`Conc(copies/µL)` * as.vector(unlist(multi_pos)) / tar_pos
  
  # write concentration into table
  name <- paste0("Concentration ", paste(genes, collapse = "."),
                 " positive for target (copies/µl)")
  if(length(conc_pos) != 0){
    tab[[name]] <- conc_pos
  } else{
    tab[[name]] <- 0
  }
  
  # set NA to zero and not used Wells (Target not in multiple positive included) to NaN
  rows_w_target <- Reduce("|", lapply(genes, function(x)(grepl(x, tab$Target))))
  tab[is.na(tab[, name]), name] = 0
  tab[!(rows_w_target), name] = NaN
  
  # compute mean concentration by target and write into table
  tab <- compute_groupwise_mean(tab, c("Target"),
                                paste0("Mean ", name),
                                name)
  
  # compute mean by (target and) well and write into table
  name2 <- paste0("Intact concentration ", paste(genes, collapse = "."), " (copies/µl)")
  tab <- compute_groupwise_mean(tab, c(), name2,
                                paste0("Mean ", name))
  
  # compute intact provirus per million cells (and shearing corrected) and write into table
  name3 <- paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))
  tab[[name3]] <- tar_mio_factor * 10^6 * tab[[name2]] / tab$`Mean concentration RPP30 + RPP30Shear (copies/µL)`
  tab[!(rows_w_target), c(name2, name3)] = NA
  tab[[paste0(name3, ", corrected for shearing")]] <- tab[[name3]] / tab$`Mean unsheared`
  
  # return extended table
  return(tab)
}

