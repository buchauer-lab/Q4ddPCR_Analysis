# file containing functions that read files (keep as short as possible!)

# read excel file
read_xlsx <- function(filename){
  library(readxl)
  
  # check if file exists
  if(!file.exists(filename)){
    stop(paste0("Specified xlsx file ", filename, " does not exist. Please check if the path and name 
         are correct and that there are no spelling mistakes."))
  }
  
  # read xlsx
  in_xlsx <- read_xlsx(filename, sheet = 1)
  
  # return
  return(in_xlsx)
}


# read csv file
read_csv <- function(filename, csv_skip, csv_nrows){
  # check if file exists
  if(!file.exists(filename)){
    stop(paste0("Specified csv file ", filename, " does not exist. Please check if the path and name 
         are correct and that there are no spelling mistakes."))
  }
  
  # read csv
  in_csv <- read.csv(csv_file, skip = csv_skip, row.names = NULL, 
                     nrows = csv_nrows)
  
  # remove shift in column names and resulting empty column at end
  colnames(in_csv) <- colnames(in_csv)[2:ncol(in_csv)]
  in_csv <- in_csv[,1:(ncol(in_csv) -1)]
  
  # return
  return(in_csv)
}


# remove channels from files (optional)
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


# create QC table
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


# remove zero channel (optional)
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
                                !(dtQC$`Sample description 1` %in% c("h2o", "H2o", "h2O", "H2O", "water", "Water", "Wasser", "wasser"))])
  
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