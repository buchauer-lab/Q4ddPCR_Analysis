# Analysis of Q4ddPCR
This is the analysis package for the data files created with Q4ddPCR.
Thus, it depends on the the input format, including certain column names (see vignette).


The package requires R verion 4.0 or higher. It was tested for R version 4.1.1 and the latest (4.5.1), on the latest stable versions (at time of creation) of Ubuntu (24.04), macOS (15), and Windows (Windows Server 2022).

The package, called "MultiplexPCRAnalyser", can be installed via

 ```devtools::install_github("buchauer-lab/Q4ddPCR_Analysis", build_vignettes = TRUE)```

Please follow the vignette for instructions ```vignette("vignette", package="MultiplexPCRAnalyser")```. The respective data and parameter file can be found at [https://github.com/Gaebler-Lab/Q4ddPCR/tree/Q4ddPCR_v1.1.0]. For version 1.0.0, please use set_parameter.R, and CSV-file.csv & Excel-file.xlsx as data input files. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15791355.svg)](https://doi.org/10.5281/zenodo.15791355)
