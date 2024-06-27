# Function to extract entries containing exactly n plus symbols
extract_exact_n_plus <- function(vector, n) {
  # Create the regular expression dynamically
  pattern <- paste0("^[^+]*", paste(rep("\\+[^+]*", n), collapse = ""), "$")
  
  # Use grepl with the generated pattern
  result <- vector[grepl(pattern, vector)]
  
  return(result)
}


get_multi_pos <- function(tab, genes, num_target = NULL){
  if(is.null(genes)){
    if(is.null(num_target)){
      stop("Need to specify genes or num target to compute multi positives.")
    }
    pattern = paste(c(rep(".*\\+", num_target), ".*"))
  } else{
    pattern <- paste(paste0("(?=.*", genes, "\\+)"), collapse = "")
    pattern <- extract_exact_n_plus(pattern, length(genes))
  }
  
  multi_pos_names <- pos_names <- grep(pattern, names(tab), value = TRUE, perl = TRUE)
  tar_pos <- unlist(apply(tab, 1, sum_target_positive_values))
  multi_pos <- apply(tab[, quad_pos_names], 1, sum)
  
  
  # TODO: change names to be variable
  name <- paste0("Concentration ", paste(genes, collapse = "."),
                 " positive for target (copies/µl)")
  tab[[name]] <- tab$`Conc(copies/µL)` * as.vector(unlist(quad_pos)) / tar_pos
  
  # do not compute for RU5 wells
  tab[!(tab$Target %in% genes), name] = NA
  
  tab <- compute_groupwise_mean(tab, c("Target"),
                                paste0("Mean ", name),
                                name)
  
  name2 <- paste0("Intact concentration ", paste(genes, collapse = "."), " (copies/µl)")
  tab <- compute_groupwise_mean(tab, c(""), name2, paste0("Mean ", name))
  
  name3 <- paste0("intact provirus/Mio cells ", paste(genes, collapse = "."))
  tab[[name3]] <- 2 * 10^6 * tab[[name2]] / tab$`Mean concentration RPP30 + RPP30Shear (copies/µL)`
  
  tab[[paste0(name3, ", corrected for shearing")]] <- tab[[name3]] / tab$`Mean unsheared`
  
  return(tab)
}

