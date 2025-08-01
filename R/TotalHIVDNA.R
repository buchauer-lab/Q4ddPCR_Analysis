#' Compute total HIV content
#'
#' Compute the total HIV content as proposed by Rachel Scheck.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @return Data frame df updated with new columns.
compute_total_HIV <- function(df) {
  # find highest Psi Env multiple (quadrupel > EnvGagPsi > EnvPsiPol > EnvPsi > EnvGagPol > EnvGag)
  # find name of column first
  name <- NULL
  if (any(grep("\\+.*\\+.*\\+.*\\+", names(df)))) {
    name <- grep("^intact provirus/Mio cells(?=.*Pol)(?=.*Gag)(?=.*Psi)(?=.*Env)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  } else if (any(grep("(?=.*Env\\+)(?=.*Gag\\+)(?=.*Psi\\+)", names(df), perl = TRUE))) {
    name <- grep("^intact provirus/Mio cells(?=.*Gag)(?=.*Psi)(?=.*Env)(?!.*Pol)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  } else if (any(grep("(?=.*Env\\+)(?=.*Psi\\+)(?=.*Pol\\+)", names(df), perl = TRUE))) {
    name <- grep("^intact provirus/Mio cells(?=.*Pol)(?=.*Psi)(?=.*Env)(?!.*Gag)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  } else if (any(grep("(?=.*Env\\+)(?=.*Psi\\+)", names(df), perl = TRUE))) {
    name <- grep("^intact provirus/Mio cells(?=.*Psi)(?=.*Env)(?!.*Pol)(?!.*Gag)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  } else if (any(grep("(?=.*Env\\+)(?=.*Gag\\+)(?=.*Pol\\+)", names(df), perl = TRUE))) {
    name <- grep("^intact provirus/Mio cells(?=.*Pol)(?=.*Gag)(?=.*Env)(?!.*Psi)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  } else if (any(grep("(?=.*Env\\+)(?=.*Gag\\+)", names(df), perl = TRUE))) {
    name <- grep("^intact provirus/Mio cells(?=.*Gag)(?=.*Env)(?!.*Psi)(?!.*Pol)(?=.*shearing)",
      names(df),
      value = TRUE, perl = TRUE
    )
  }
  # then extract value
  if (is.null(name)) {
    quad <- 0
  } else if(length(name) > 1){
    stop("More than 1 column for Total HIV calculation selected.")
  } else{
    quad <- na.omit(unique(df[[name]]))
    # check if length of value is correct
    if(length(quad) == 0){
      warning("Quadruple positives were all NA and are set to 0.")
      quad <- 0
    } else if (length(quad) > 1) {
      stop("Bug: length of quad must be one. Please contact maintainer.")
    }
  }
  substrings <- c("Env", "Gag", "Pol", "Psi")
  # get columns with 3 substrings
  columns_with_3_substrings <- grep("^intact(?!.*shearing)", names(df)[sapply(names(df), function(name) {
    # Count how many substrings are present in the column name
    sum(sapply(substrings, grepl, name)) == 3
  })], value = TRUE, perl = T)
  # compute triplet sum
  trip <- sum(apply(df[, columns_with_3_substrings], 2, unique), na.rm = T)

  # get columns with 2 substrings
  columns_with_2_substrings <- grep("^intact(?!.*shearing)", names(df)[sapply(names(df), function(name) {
    # Count how many substrings are present in the column name
    sum(sapply(substrings, grepl, name)) == 2
  })], value = TRUE, perl = T)
  # compute triplet sum
  doub <- sum(apply(df[, columns_with_2_substrings], 2, unique), na.rm = T)

  # compute for singlets
  sing <- sum(df[["Mean Target/Mio cells"]]) / 4

  # Compute total HIV DNA / Mio cells
  df[["total HIV DNA/Mio cells"]] <- sing - doub + trip - quad
  return(df)
}

#' Compute total HIV content
#'
#' Estimate the total HIV content by Env Psi content.
#' @param df The data frame for whose column the total HIV content shall be computed
#' @import dplyr
#' @return Data frame df updated with new column.
compute_total_HIV_envPsi <- function(df) {
  # check for column names
  if(!("Mean Target/Mio cells" %in% names(df))){
    stop("No Mean Target/Mio cells detected")
  }
  
  # extract ENV, Psi, und EnvPsi concentration from df (per Mio cells)
  env <- unique(df[grepl("Env", df$Target), "Mean Target/Mio cells"])
  psi <- unique(df[grepl("Psi", df$Target), "Mean Target/Mio cells"])
  if(!any(c("intact provirus/Mio cells Psi.Env, corrected for shearing",
            "intact provirus/Mio cells Env.Psi, corrected for shearing") 
          %in% names(df))){
    warning("No Env+Psi+ detected. Will use 0.")
    envpsi <- 0
  } else{
    envpsi <- unique(na.omit(df$`intact provirus/Mio cells Psi.Env, corrected for shearing`))
    if(is.null(envpsi)){
      envpsi <- unique(na.omit(df$`intact provirus/Mio cells Env.Psi, corrected for shearing`))
    }
  }
  df[["total HIV DNA/Mio cells (Env.Psi)"]] <- unlist(unname(env + psi - envpsi))
  return(df)
}
