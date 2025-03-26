# file containing functions that read files (keep as short as possible!)

#' Read xlsx file
#'
#' Read the excel file (first sheet) of the input files
#' @param filename Path to an existing xlsx file.
#' @return Data frame with the contents of the xlsx file. NA values in
#' concentration column are filled in with zeros.
read_xlsx <- function(filename) {
  # check if file exists
  if (!file.exists(filename)) {
    stop(paste0("Specified xlsx file ", filename, " does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  }

  # read xlsx
  in_xlsx <- openxlsx::read.xlsx(filename, sheet = 1, sep.names = " ")

  # adapt concentration column
  suppressWarnings(in_xlsx$`Conc(copies/µL)` <- as.numeric(in_xlsx$`Conc(copies/µL)`))
  in_xlsx$`Conc(copies/µL)`[is.na(in_xlsx$`Conc(copies/µL)`)] <- 0

  # return
  return(in_xlsx)
}


#' Read csv file
#'
#' Read the csv file of the input files within the specified rows.
#' @param filename Path to an existing csv file.
#' @param csv_skip Number of lines in csv file that should not be read in, since
#' they are before the actual table. Otherwise reading the file will return an
#' error due to column number mismatch in different rows.
#' @param csv_end Last row number in csv file with relevant data. Otherwise the
#' following table (in the expected input) will cause an error due to column
#' number mismatch in different rows.
#' @return Data frame with the relevant data.
read_csv <- function(filename, csv_skip) {
  # check if file exists
  if (!file.exists(filename)) {
    stop(paste0("Specified csv file ", filename, " does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  }

  # read csv
  in_csv <- data.table::fread(filename, skip = csv_skip, fill = TRUE)

  # remove second table
  if (is.null(in_csv$Well)) {
    stop("csv file must contain column Well.")
  }

  # ensure that csv structure is as expected (two tables in one file)
  if (!(is.na(which(in_csv$Well == "")[1] - 1))) {
    in_csv <- in_csv[1:which(in_csv$Well == "")[1] - 1, ]
  } else {
    warning("csv deviates from expected format. Trying to continue.")
  }


  # remove empty column due to commas at end of file
  in_csv <- in_csv[, 1:(ncol(in_csv) - 1)]

  # check if correct number of rows was read
  if (any(in_csv$Well == "Well")) {
    l <- 1 + length(in_csv$Well) - which(in_csv$Well == "Well")
    stop(paste0("Too many rows from csv are read. Potentially ", l, " too many."))
  }

  # return
  return(in_csv)
}


#' Remove channels manually (optional)
#'
#' Remove channels if they are not clean or irrelevant.
#' @param df Data frame from which the channels are to be removed. Usually the
#' output of reading the csv and the xlsx file each.
#' @param channels Channels to remove (Well ID)
#' @return Data frame without the specified channels
rm_channels <- function(df, channels) {
  # check if df has Well column
  if (is.null(df$Well)) {
    stop("No Wells specified in data, thus cannot be removed. Please check data.")
  }
  # check if all channels are present
  if (any(!(channels %in% df$Well))) {
    stop("At least one of the Wells that should be removed does not occur in
         the Well column of the csv file. Please only remove Wells that do occur
         in the file.")
  }
  # remove
  df <- df[!(df$Well %in% channels), ]

  # return
  return(df)
}


#' Create QC table
#'
#' Create the dtQC table, which contains the relevant content from the xlsx file.
#' @param xlsx_df The data frame created when reading the xlsx file.
#' @param cols_of_int Columns that should be contained in dtQC. If not specified,
#' a standard set of columns is kept. (optional)
#' @param acc_drop_factor Numeric factor by which the threshold is computed from
#' the accepted droplets. (default=3)
#' @return The dtQC dataframe with the relevant columns.
create_dtQC <- function(xlsx_df, cols_of_int = NULL, acc_drop_factor = 3) {
  if (is.null(cols_of_int)) {
    cols_of_int <- c(
      "Well", "Sample description 1", "DyeName(s)", "Target",
      "Conc(copies/µL)", "Accepted Droplets", "Positives",
      "Negatives", grep("Ch", names(xlsx_df), value = TRUE)
    )
  }

  # Ensure all required columns exist
  if (!all(cols_of_int %in% names(xlsx_df))) {
    stop("Some required columns are missing from the data frame.")
  }

  dtQC <- xlsx_df[, cols_of_int]
  dtQC$Threshold <- dtQC$`Accepted Droplets` / acc_drop_factor
  dtQC$`Total positives` <- apply(dtQC[, grep("\\+", names(dtQC), value = TRUE)], 1, sum)
  return(dtQC)
}


#' Remove zero channel
#'
#' Remove wells for which at least one channel is defect, i.e., has
#' concentration 0. Wells with H2O probes are not removed.
#' @param dtQC Data frame created with create_dtQC.
#' @param in_csv Data frame created with read_csv.
#' @return List of data frames dtQC and in_csv without the defect channels.
rm_zero_channel <- function(dtQC, in_csv) {
  # check if expected columns exist
  if (is.null(dtQC$`Conc(copies/µL)`)) {
    stop("No concentration column 'Conc(copies/µl)' found in dtQC. Please run create_dtQC first.")
  }
  if (is.null(dtQC$`Sample description 1`)) {
    stop("No 'Sample description 1' column found in dtQC. Please run create_dtQC first.")
  }

  # get zero channels without H20 channels
  zero_ch <- unique(dtQC$Well[(dtQC$`Conc(copies/µL)` == 0) &
    !(dtQC$`Sample description 1` %in% c("h2o", "H2o", "h2O", "H2O"))])

  # make warning if there are zero channels
  if (length(zero_ch) > 0) {
    warning(paste0(
      "These wells will be removed, as they have at least one channel with concentration value 0: ",
      paste(zero_ch, collapse = ", ")
    ))
  }

  # remove channels from data objects
  dtQC <- dtQC[!(dtQC$Well %in% zero_ch), ]
  in_csv <- in_csv[!(in_csv$Well %in% zero_ch), ]

  # return
  return(list(dtQC, in_csv))
}

#' Read files
#'
#' Read the input xlsx and csv file to create the initial tables for the analysis.
#' @param xlsx_file Path to an existing xlsx_file.
#' @param csv_file Path to an existing csv_file.
#' @param csv_skip Number of lines in csv file that should not be read in, since
#' they are before the actual table. Otherwise reading the file will return an
#' error due to column number mismatch in different rows.
#' @param remove_channel List of channels that are to be removed. (optional)
#' @param rm_zero_channels Boolean value to indicate whether defect channels,
#' i.e., channels with zero concentration are to be removed. (optional)
#' @export
read_files <- function(xlsx_file,
                       csv_file,
                       csv_skip,
                       remove_channel = NULL,
                       rm_zero_channels = FALSE) {
  # read files
  in_csv <- read_csv(csv_file, csv_skip)
  in_xlsx <- read_xlsx(xlsx_file)

  # remove specified channels
  if (!(is.null(remove_channel))) {
    in_csv <- rm_channels(in_csv, remove_channel)
    in_xlsx <- rm_channels(in_xlsx, remove_channel)
  }


  # create sheet Data_table Quality Check (dtQC)
  dtQC <- create_dtQC(in_xlsx)

  # remove wells that have concentration zero for at least one well if specified
  if (rm_zero_channels) {
    out.list <- rm_zero_channel(dtQC, in_csv)
    return(out.list)
  } else {
    return(list(in_csv, dtQC))
  }
}

#' Read multiple files
#'
#' Read a list of input xlsx and csv files to create the initial tables for the analysis.
#' @param xlsx_files Vector with paths to existing xlsx_files.
#' @param csv_files Vector with paths to existing csv_files.
#' @param csv_skip Number of lines in csv file that should not be read in, since
#' they are before the actual table. Otherwise reading the file will return an
#' error due to column number mismatch in different rows. Same in every file.
#' @param remove_channel List of lists of channels that are to be removed. Will
#' be processed in the same order as the csv and xlsx files. (optional)
#' @param rm_zero_channels Boolean value to indicate whether defect channels,
#' i.e., channels with zero concentration are to be removed. (optional)
#' @export
read_multiple_files <- function(xlsx_files,
                                csv_files,
                                csv_skip = 4,
                                remove_channel = NULL,
                                rm_zero_channels = FALSE) {
  # check that files match
  if (length(xlsx_files) != length(csv_files)) {
    stop("There must be as many xlsx files as there are csv files.")
  }

  # lists to save outputs
  all_csv <- list()
  all_dtQC <- list()

  # save names of data frames
  csv_names <- list()
  dtQC_names <- list()

  # read each pair of files
  for (i in 1:length(xlsx_files)) {
    if (is.null(remove_channel)) {
      tmp <- read_files(xlsx_files[i], csv_files[i], csv_skip, rm_zero_channels = rm_zero_channels)
    } else {
      tmp <- read_files(xlsx_files[i], csv_files[i], csv_skip,
        remove_channel = remove_channel[i],
        rm_zero_channels = rm_zero_channels
      )
    }

    # get output
    in_csv <- tmp[[1]]
    dtQC <- tmp[[2]]

    # Make wells unique by adding plate number
    in_csv$Well <- paste0(in_csv$Well, "_", i)
    dtQC$Well <- paste0(dtQC$Well, "_", i)

    # save in lists
    all_csv[[i]] <- in_csv
    all_dtQC[[i]] <- dtQC

    csv_names[[i]] <- colnames(in_csv)
    dtQC_names[[i]] <- colnames(dtQC)
  }
  # check that all lists have the same names
  if (!(all(sapply(csv_names, function(x) identical(x, csv_names[[1]]))))) {
    stop("CSV files have different column names.")
  }
  if (!(all(sapply(dtQC_names, function(x) identical(x, dtQC_names[[1]]))))) {
    stop("CSV files have different column names.")
  }
  # combine lists
  final_csv <- do.call(rbind, all_csv)
  final_dtQC <- do.call(rbind, all_dtQC)

  return(list(final_csv, final_dtQC))
}
