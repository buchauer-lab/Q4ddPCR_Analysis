# branch dev
# clear workspace
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source("ExtractInfos.R")

# set arguments

csv_file <- "Data5/14.05.2024_ClusterData 1.csv"
csv_skip <- 4 # number of rows before the table starts
csv_nrows <- 389 - (csv_skip + 1) # number of rows in the table
xlsx_file <- "Data5/14.05.2024_DataSheet 1.xlsx"
output_file <- "Data5/Output.xlsx"

ch_dye <- c("Ch1" = "FAM",
            "Ch2" = "VIC",
            "Ch3" = "Cy5",
            "Ch5" = "ROX",
            "Ch6" = "ATTO590")


custom_dilution_factor <- FALSE
dilution_factor <- c("gDNA 8E5 0,3 ng" = 100,
                     "gDNA 8E5 3 ng" = 100,
                     "gDNA 8E5 30 ng" = 100)

remove_channel <- c()

# ================== execute functions =========================================
# read files
information <- read_files(xlsx_file, csv_file, csv_skip, csv_nrows, remove_channel)
in_csv <- information[[1]]
dtQC <- information[[2]]

# create standard table and values
standards <- create_household_table(dtQC, dilution_factor, custom_dilution_factor)
tab1 <- standards[[1]]
mean_conc_household <- standards[[2]]
mean_unsheared <- standards[[3]]
mean_cells_per_reac <- standards[[4]]
Shear_ind <- standards[[5]]
rpp30_names <- standards[[6]]
df_tw <- standards[[7]]

# define groups to be analyzed together
split_into_groups <- define_groups(df_tw, rpp30_names)
mat_groups <- split_into_groups[[1]]
groups <- split_into_groups[[2]]

# create output tables
output <- create_tables(mat_groups, in_csv, groups, ch_dye)
output_tables <- output[[1]]
conf_mats <- output[[2]]

# write to xlsx file
write_output_file(output_tables, conf_mats, tab1, output_file) 











