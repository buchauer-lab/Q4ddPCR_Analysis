# file containing functions that read files (keep as short as possible!)

#' Read xlsx file
#' 
#' Read the excel file (first sheet) of the input files
#' @param filename Path to an existing xlsx file.
#' @return Data frame with the contents of the xlsx file. NA valuess in
#' concentration column are filled in with zeros.
#' @export
read_xlsx <- function(filename){
  library(readxl)
  
  # check if file exists
  if(!file.exists(filename)){
    stop(paste0("Specified xlsx file ", filename, " does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  }
  
  # read xlsx
  in_xlsx <- readxl::read_xlsx(filename, sheet = 1)
  
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
#' @export
read_csv <- function(filename, csv_skip, csv_end){
  # check if file exists
  if(!file.exists(filename)){
    stop(paste0("Specified csv file ", filename, " does not exist. Please check if the path and name are correct and that there are no spelling mistakes."))
  }
  
  # compute csv nrows
  csv_nrows <- csv_end - (csv_skip + 1)
  
  # read csv
  in_csv <- read.csv(filename, skip = csv_skip, row.names = NULL, 
                     nrows = csv_nrows)
  
  # remove shift in column names and resulting empty column at end
  colnames(in_csv) <- colnames(in_csv)[2:ncol(in_csv)]
  in_csv <- in_csv[,1:(ncol(in_csv) -1)]
  
  # check if correct number of rows was read
  if(any(in_csv$Well == "Well")){
    l = 1 + length(in_csv$Well) - which(in_csv$Well == "Well")
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
#' @export
rm_channels <- function(df, channels){
  
  # check if df has Well column
  if(is.null(df$Well)){
    stop("No Wells specified in data, thus cannot be removed. Please check data.")
  }
  # check if all channels are present
  if(any(!(channels %in% df$Well))){
    stop("At least one of the Wells that should be removed does not occur in
         the Well column of the csv file. Please only remove Wells that do occur
         in the file.")
  }
  # remove
  df <- df[!(df$Well %in% channel),]
  
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
#' @export
create_dtQC <- function(xlsx_df, cols_of_int = NULL, acc_drop_factor = 3){
  
  # add cols_of_int if necessary
  if(is.null(cols_of_int)){
  cols_of_int <- c("Well", "Sample description 1", "DyeName(s)", "Target", 
                   "Conc(copies/µL)", "Accepted Droplets", "Positives", 
                   "Negatives", grep("Ch", names(in_xlsx), value = T))
  }
  
  # create sheet Data_table Quality Check (dtQC)
  dtQC <- xlsx_file[,cols_of_int]
  
  # add threshold and number of total positives
  dtQC$Threshold <- dtQC$`Accepted Droplets`/acc_drop_factor
  dtQC$`Total positives` <- apply(dtQC[,grep("\\+", names(dtQC), value=T)], 1, sum)
  
  return(dtQC)
}


#' Remove zero channel 
#' 
#' Remove wells for which at least one channel is defect, i.e., has
#' concentration 0. Wells with H2O probes are not removed.
#' @param dtQC Data frame created with create_dtQC.
#' @param in_csv Data frame created with read_csv.
#' @return List of data frames dtQC and in_csv without the defect channels.
#' @export
rm_zero_channel <- function(dtQC, in_csv){
  
  # check if expected columns exist
  if(is.null(dtQC$`Conc(copies/µL)`)){
    stop("No concentration column 'Conc(copies/µl)' found in dtQC. Please run create_dtQC first.")
  }
  if(is.null(dtQC$`Sample description 1`)){
    stop("No 'Sample description 1' column found in dtQC. Please run create_dtQC first.")
  }
  
  # get zero channels without H20 channels
  zero_ch <- unique(dtQC$Well[(dtQC$`Conc(copies/µL)` == 0) & 
                                !(dtQC$`Sample description 1` %in% c("h2o", "H2o", "h2O", "H2O"))])
  
  # make warning if there are zero channels
  if(length(zero_ch) > 0){
    warning(paste0("These wells will be removed, as they have at least one channel with concentration value 0: ",
                   paste(zero_ch, collapse = ", ")))
  }
  
  # remove channels from data objects
  dtQC <- dtQC[!(dtQC$Well %in% zero_ch),]
  in_csv <- in_csv[!(in_csv$Well %in% zero_ch),]
  
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
#' @param csv_end Last row number in csv file with relevant data. Otherwise the
#' following table (in the expected input) will cause an error due to column
#' number mismatch in different rows.
#' @param remove_channel List of channels that are to be removed. (optional)
#' @param rm_zero_channels Boolean value to indicate whether defect channels, 
#' i.e., channels with zero concentration are to be removed. (optional)
#' @export
read_files <- function(xlsx_file,
                       csv_file,
                       csv_skip,
                       csv_nrows,
                       remove_channel = NULL,
                       rm_zero_channels = FALSE){
  
  # read files
  in_csv <- read_csv(csv_file, csv_skip, csv_nrows)
  in_xlsx <- read_xlsc(xlsx_file)
  
  # remove specified channels
  if(!(is.null(remove_channel))){
    in_csv <- rm_channels(in_csv, remove_channel)
    in_xlsx <- rm_channels(in_xlsx, remove_channel)
  }

  
  # create sheet Data_table Quality Check (dtQC)
  dtQC <- create_dtQC(in_xlsx)
  
  # remove wells that have concentration zero for at least one well if specified
  if(rm_zero_channels){
    out.list <- rm_zero_channels(dtQC, in_csv)
    return(out.list)
  } else{
    return(list(in_csv, dtQC))
  }
}