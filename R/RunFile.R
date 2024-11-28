
# clear workspace
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source("ExtractInfos.R")

# set arguments


csv_file <- "../../data/Data16/20241015_GTG_Q4_A5321_Plate2_20241015_213719_533_ClusterData.csv"
csv_skip <- 4 # number of rows before the table starts
csv_nrows <- 432 - (csv_skip + 1) # number of rows in the table
xlsx_file <- "../../data/Data16/20241015_GTG_Q4_A5321_Plate2_20241015_213719_533.xlsx"
output_file <- "../../data/Data16/Output.xlsx"

ch_dye <- c("Ch1" = "FAM",
            "Ch2" = "VIC",
            "Ch3" = "Cy5",
            "Ch5" = "ROX",
            "Ch6" = "ATTO590")

custom_dilution_factor <- FALSE

dilution_factor <- c("112213 TP1" = 100,
                     "112213 tp2" = 100,
                     "214356" = 100,
                     "232656" = 100,
                     "581105" = 100)
#dilution_factor <- c("gDNA 5104" = 100, data3
#                     "gDNA 8E5" = 1)

remove_channel <- c("A10", "B10", "C10", "D10", "E10", "F10", "G10", "H10",
                    "A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11",
                    "A09", "B09", "C09", "D09", "B07", "C07", "D07")
rm_zero_channel_wells <- F # remove wells that have concentration 0 for at least one channel
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

# create standard table and values
standards <- create_household_table(dtQC, dilution_factor, custom_dilution_factor, threshold, mean_copies_factor, mean_cells_per_reac_factor)
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
output <- create_tables(mat_groups, in_csv, groups, ch_dye, multi_positives, threshold, tar_mio_factor)
output_tables <- output[[1]]
conf_mats <- output[[2]]
h2o_tables <- output[[3]]

# write to xlsx file
write_output_file(output_tables, conf_mats, tab1, output_file, h2o_tables, multi_positives) 











 

