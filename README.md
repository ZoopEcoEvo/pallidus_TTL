# Title	

Author<sup>1</sup>, Author<sup>2</sup>, Author<sup>3</sup> 	

1. Affiliation	
2. Affiliation	
3. Affiliation	

Project description. Includes basic info about collection locations and times. 


<a href="DOI LINK GOES HERE"><img src="https://img.shields.io/badge/Publication%20Here-Here?style=for-the-badge&labelColor=grey&color=577D58"/></a>    <a href="PREPRINT LINK GOES HERE"><img src="https://img.shields.io/badge/Preprint%20Here-Here?style=for-the-badge&labelColor=grey&color=6F92B8"/></a>    <a href="DATA LINK GOES HERE"><img src="https://img.shields.io/badge/Dataset%20Here-Here?style=for-the-badge&labelColor=grey&color=EAB221"/></a>

## Directory Structure 
The root directory contains the README and Licensing files, along with a .Rproj file and four sub-directories: Data, Manuscript, Output, and Scripts.  

-   `Data/` contains the raw data used in this analysis.  

-   `Manuscript/` contains the R Markdown file, templates, and bibliography used to produce the manuscript draft. 

-   `Output/` contains the various products of the project (processed data, figures, knit reports, and a PDF copy of the manuscript). Note, the `Reports/` directory contains the R Markdown file used to generate the figures used in the manuscript.  

-   `Scripts/` contains two R scripts. 
    -   `01_Data_analysis.R` is used to process the raw data. The primary components of this analysis are ____. 
    -   `02_make_report.R` is use to control the project workflow. Through this script, you can choose to run the process the data, make the figures, or knit the manuscript. This script should be used rather than running isolated fragments individually. 


## Data Description 	

The `Data/` directory contains the ____ data sets used in this study.  	

-   `data.xlsx` contains data from egg production experiments during the simulated heatwave experiments for the field-collected females.   	
    -   *col_name1* - Description	  	
    -   *col_name2*	- Description 	
    -   *col_name3* - Description		   


## Workflow

The workflow is operated via the 02_Make_report.R script in the Scripts directory. It is not recommended that you run analyses or knit documents from the files themselves as the file paths are internally set and may not be correct otherwise. At the top of this script, you are able to indicate whether:

1. The raw data should be processed to ____.  

2. The summary file (located in the Output/Reports directory) should be knit. This markdown file will generate the figures used in the manuscript, as well as an HTML and a GitHub flavored markdown document.

3. The manuscript file (located in the Manuscripts directory) should be knit. This markdown file will produce formatted PDF and word document versions of the manuscript. 


## Versioning   

R version 4.2.2 (2022-10-31)  

Platform: x86_64-apple-darwin17.0 (64-bit)  

Running under: macOS Ventura 13.3.1 
  
**Attached base packages:** stats, graphics, grDevices, utils, datasets, methods, base     

**Other attached packages:** forcats_0.5.2, stringr_1.5.0, dplyr_1.0.10, purrr_1.0.0, readr_2.1.3, tidyr_1.2.1, tibble_3.1.8, tidyverse_1.3.2, ggrepel_0.9.2, ggpubr_0.5.0, ggplot2_3.4.0, car_3.1-1, carData_3.0-5, dabestr_0.3.0, magrittr_2.0.3, minpack.lm_1.2-2, MASS_7.3-58.1, boot_1.3-28, broom_1.0.2, nls.multstart_1.2.0, rTPC_1.0.2, readxl_1.4.1, knitr_1.41, rmarkdown_2.20, Rcpp_1.0.10, lubridate_1.9.0, assertthat_0.2.1, digest_0.6.31, utf8_1.2.2, R6_2.5.1, cellranger_1.1.0, backports_1.4.1, reprex_2.0.2, evaluate_0.19, httr_1.4.4, pillar_1.8.1, rlang_1.0.6, googlesheets4_1.0.1, rstudioapi_0.14, googledrive_2.0.0, munsell_0.5.0, compiler_4.2.2, modelr_0.1.10, xfun_0.36, pkgconfig_2.0.3, htmltools_0.5.4, tidyselect_1.2.0, fansi_1.0.3, crayon_1.5.2, tzdb_0.3.0, dbplyr_2.2.1, withr_2.5.0, grid_4.2.2, jsonlite_1.8.4, gtable_0.3.1, lifecycle_1.0.3, DBI_1.1.3, scales_1.2.1, cli_3.5.0, stringi_1.7.8, ggsignif_0.6.4, fs_1.5.2, xml2_1.3.3, ellipsis_0.3.2, generics_0.1.3, vctrs_0.5.1, tools_4.2.2, glue_1.6.2, hms_1.1.2, rsconnect_0.8.28, abind_1.4-5, fastmap_1.1.0, timechange_0.1.1, colorspace_2.0-3, gargle_1.2.1, rvest_1.0.3, rstatix_0.7.1, haven_2.5.1


## Funding

This study was funded by grants from ____.
