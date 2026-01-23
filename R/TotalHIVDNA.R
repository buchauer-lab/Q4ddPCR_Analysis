#' Get Quad-pos. column
#'
#' get column to use as shear-corrected quadrupel positives with the order 
#' (quadrupel > EnvGagPsi > EnvPsiPol > EnvPsi > EnvGagPol > EnvGag)
#' @param df data table
#' @import stats
#' @return quad positive value
get_quad_pos <- function(df){
  # find highest Psi Env multiple (quadrupel > EnvGagPsi > EnvPsiPol > EnvPsi > EnvGagPol > EnvGag)
  # find name of column first
  quad_order <- c("^intact provirus/Mio cells(?=.*Pol)(?=.*Gag)(?=.*Psi)(?=.*Env)(?=.*shearing)", # quad positive
                  "^intact provirus/Mio cells(?=.*Gag)(?=.*Psi)(?=.*Env)(?!.*Pol)(?=.*shearing)", # Gag+Psi+Env+
                  "^intact provirus/Mio cells(?=.*Pol)(?=.*Psi)(?=.*Env)(?!.*Gag)(?=.*shearing)", # Pol+Psi+Env+
                  "^intact provirus/Mio cells(?=.*Psi)(?=.*Env)(?!.*Pol)(?!.*Gag)(?=.*shearing)", # Psi+Env+
                  "^intact provirus/Mio cells(?=.*Pol)(?=.*Gag)(?=.*Env)(?!.*Psi)(?=.*shearing)", # Pol+Gag+Env+
                  "^intact provirus/Mio cells(?=.*Gag)(?=.*Env)(?!.*Psi)(?!.*Pol)(?=.*shearing)"  # Gag+Env+
                  )
  quad <- 0
  for (comb in quad_order) {
    name <- grep(comb,
                 names(df),
                 value = TRUE, perl = TRUE)
    print(name)
    if(length(name) == 1){
      quad <- na.omit(unique(df[[name]]))
      if((length(quad) != 0) && (quad != 0))
        break # found quad, leave for loop
    }
  }
  print(quad)
  return(quad)
}

#' Get multiplet column
#'
#' get columns to use for multiplet (triplet or doublet) computation
#' @param df data table
#' @param multiplet number of targets to look for (2 for doublets, 3 for multiplets)
#' @param substrings list of targets available
#' @return sum of multiplets
get_intact_multiplet_conc <- function(df, multiplet, substrings){
  # get columns with n substrings
  columns_with_n_substrings <- grep("^intact(?!.*shearing)", names(df)[sapply(names(df), function(name) {
    # Count how many substrings are present in the column name
    sum(sapply(substrings, grepl, name)) == multiplet
  })], value = TRUE, perl = T)
  # compute triplet sum
  return(sum(apply(df[, columns_with_n_substrings], 2, unique), na.rm = T))
}

#' Compute total HIV content
#'
#' Compute the total HIV content. If quadruple positive values are not found in the
#' data, the following order of multiplets will be used until a substitute is found:
#' (quadrupel > EnvGagPsi > EnvPsiPol > EnvPsi > EnvGagPol > EnvGag). For more information
#' check the associated publication.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @return Data frame df updated with new columns.
#' @export
compute_total_HIV <- function(df) {
  # get value for quadruple positives
  quad <- get_quad_pos(df)

  substrings <- c("Env", "Gag", "Pol", "Psi") # TODO: get from data

  # compute triplet sum
  trip <- get_intact_multiplet_conc(df,  3, substrings)
  # compute doublet sum
  doub <- get_intact_multiplet_conc(df,  2, substrings)
  # compute for singlets
  sing <- sum(df[["Mean Target/Mio cells"]])/length(unique(df$Well)) # as mean is taken for each target over all Wells: divide by number of Wells
  # Compute total HIV DNA / Mio cells
  return(sing - doub + trip - quad)
}

#' Compute total HIV content
#'
#' Estimate the total HIV content by Env Psi content.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @import stats
#' @return Data frame df updated with new column.
#' @export
compute_total_HIV_envPsi <- function(df) {
  # check for column names
  if(!("Mean Target/Mio cells" %in% names(df))){
    stop("No Mean Target/Mio cells detected")
  }
  
  # extract ENV, Psi, und EnvPsi concentration from df (per Mio cells)
  env <- unique(df[grepl("Env", df$Target), "Mean Target/Mio cells"])
  psi <- unique(df[grepl("Psi", df$Target), "Mean Target/Mio cells"])
  
  # get column name with Env+Psi+
  name <- grep("^intact provirus/Mio cells(?=.*Psi)(?=.*Env)(?!.*Pol)(?!.*Gag)(?=.*shearing)",
               names(df),
               value = TRUE, perl = TRUE)
  if(length(name) == 0){
    warning("No Env+Psi+ detected. Will use 0.")
    envpsi <- 0
  } else{
    envpsi <- unique(na.omit(df[[name]]))
  }
  
  # return "total HIV DNA/Mio cells (Env.Psi)"
  return(unlist(unname(env + psi - envpsi)))
}
