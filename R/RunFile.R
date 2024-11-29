
# clear workspace
source("R/ExtractInfos.R")
source("R/CreateHouseholdTable.R")
source("R/Readers.R")
setwd("R/")
setwd("..")

# set arguments
csv_file <- "../../data/Data14/20240807_GTG_.csv"
csv_skip <- 4 # number of rows before the table starts
csv_nrows <- 261 - (csv_skip + 1) # number of rows in the table
xlsx_file <- "../../data/Data14/20240807_GTG_.xlsx"
output_file <- "../../data/Data14/Output.xlsx"

ch_dye <- c("Ch1" = "FAM",
            "Ch2" = "VIC",
            "Ch3" = "Cy5",
            "Ch5" = "ROX",
            "Ch6" = "ATTO590")

custom_dilution_factor <- TRUE

dilution_factor <- c("112213 TP1" = 100,
                     "112213 tp2" = 100,
                     "214356" = 100,
                     "232656" = 100,
                     "581105" = 100)
#dilution_factor <- c("gDNA 5104" = 100,
 #                    "gDNA 8E5" = 1)

remove_channel <- c("A04","B04","C04","D04","C01","D01","C02","D02","C06","D06","C05","D05","A07","B07","C07","D07","E07","F07","G07","H07","E06","F06","G06","H06")
rm_zero_channel_wells <- FALSE # remove wells that have concentration 0 for at least one channel
                              # will not remove H2O channels

compute_all_positives_for <- c("Psi", "Env", "Gag", "Pol")

multi_positives <- get_multipos(compute_all_positives_for)


# define minimum number of accepted droplets to continue with well
threshold <- 7500
mean_copies_factor <- 20 # number to multiply Mean concentration RPP30 + Shear with to compute mean copies/well # I guess should be named volume
mean_cells_per_reac_factor <- 2 # factor to multiply Mean copies/cell with to obtain Mean cells per reaction
tar_mio_factor <- 2 # factor to multiply Concentration with to obtain Target/Mio cells (same for multiple positives)
# ================== execute functions =========================================
# read files
information <- read_files(xlsx_file, csv_file, csv_skip, csv_nrows, remove_channel, rm_zero_channel_wells)
in_csv <- information[[1]]
dtQC <- information[[2]]

# group data according to sample description and targets
grouped_data <- define_groups(dtQC, dilution_factor)

# create standard table and values
standards <- create_household_table(dtQC, grouped_data, thresh, mean_copies_factor, mean_cells_per_reac_factor)
tab1 <- standards[[1]]
grouped_data <- standards[[2]]

# TODO: adapt below
# create output tables
#output <- create_tables(mat_groups, in_csv, groups, ch_dye, multi_positives, threshold, tar_mio_factor)
#output_tables <- output[[1]]
#conf_mats <- output[[2]]
#h2o_tables <- output[[3]]

# write to xlsx file
#write_output_file(output_tables, conf_mats, tab1, output_file, h2o_tables, multi_positives) 











 

