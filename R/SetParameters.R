csv_file <- "../data/Data21/20251110_Q4ddPCR_PS_20251110_164059_397_ClusterData.csv"
csv_skip <- 4 # number of rows before the table starts
xlsx_file <- "../data/Data21/20251110_Q4ddPCR_PS_20251110_164059_397.xlsx"
output_file <- "../data/Data21/Output.xlsx"

ch_dye <- c("Ch1" = "FAM",
            "Ch2" = "VIC",
            "Ch3" = "Cy5",
            "Ch5" = "ROX",
            "Ch6" = "ATTO590")

target_channel <- c("Target 1" = "Ch1",
                    "Target 2" = "Ch2",
                    "Target 3" = "Ch3",
                    "Target 4" = "Ch6")

custom_dilution_factor <- FALSE


dilution_factor <- c("iCARE 1" = 100,
                     "iCARE 2" = 100,
                     "iCARE 3" = 100,
                     "iCARE 4" = 100,
                     "iCARE 5" = 100,
                     "iCARE 6" = 100,
                     "iCARE 7" = 100
)


remove_channel <- c("E11","F11","G11","H11",
                    "A12", "B12","C12", "D12", "E12", "F12", "G12", "H12"
                    
)

# 
#"A10", "B10","C10", "D10", "E10", "F10", "G10", "H10",
#"A07", "B07","C07", "D07", "E07", "F07", "G07", "H07",
#"A08", "B08","C08", "D08", "E08", "F08", "G08", "H08",
#"A06", "B06","C06", "D06", "E06", "F06", "G06", "H06",
#
#"A10", "B10","C10", "D10", "E10", "F10", "G10", "H10"
#
#
#
rm_zero_channel_wells <- FALSE # remove wells that have concentration 0 for at least one channel
# will not remove H2O channels

compute_all_positives_for <- c("Psi", "Env", "Gag", "Pol")

multi_positives <- get_multipos(compute_all_positives_for)

tar_mio_factor <- c("iCARE 1" = 2,
                    "iCARE 2" = 2,
                    "iCARE 3" = 2,
                    "iCARE 4" = 2,
                    "iCARE 5" = 2,
                    "iCARE 6" = 2,
                    "iCARE 7" = 2)

# define minimum number of accepted droplets to continue with well
threshold <- 7500
mean_copies_factor <- 20 # number to multiply Mean concentration RPP30 + Shear with to compute mean copies/well # I guess should be named volume
mean_cells_per_reac_factor <- c("iCARE 1" = 2,
                                "iCARE 2" = 2,
                                "iCARE 3" = 2,
                                "iCARE 4" = 2,
                                "iCARE 5" = 2,
                                "iCARE 6" = 2,
                                "iCARE 7" = 2)
# ================== execute functions ===========

shear_name <- c("RPP30", "RPP30Shear")
water_name <- c("H20")