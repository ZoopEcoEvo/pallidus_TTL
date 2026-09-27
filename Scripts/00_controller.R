# Load in required packages
source("Scripts/rezende_functions.R")

library(rmarkdown)
library(survival)
library(broom)
library(scales)
library(tidyverse)

#Determine which scripts should be run
process_data = F #Runs data analysis 
make_report = T #Runs project summary
knit_manuscript = F #Compiles manuscript draft

############################
### Read in the RAW data ###
############################

if(process_data == T){
  source(file = "Scripts/01_data_processing.R")
}

files = list.files("Raw_data/surv_data/", pattern = "\\.csv$", full.names = TRUE)

raw_data = files |>
  map(read_csv, show_col_types = FALSE) |>
  list_rbind(names_to = "source_file")

env_temps = read.csv("Raw_data/stream_temps/22229810 2026-06-01 14_30_55 EDT.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::select("date_time" = date_time_edt, 
                "temp_c" = temperature_c) %>% 
  mutate(date_time = as_datetime(date_time, format = "%m/%d/%Y %H:%M:%S"))

env_temps_minutes = env_temps |> 
  uncount(10) |> 
  drop_na()

##################################
### Read in the PROCESSED data ###
##################################

if(make_report == T){
  render(input = "Output/Reports/report.Rmd", #Input the path to your .Rmd file here
         #output_file = "report", #Name your file here if you want it to have a different name; leave off the .html, .md, etc. - it will add the correct one automatically
         output_format = "all")
}

##################################
### Read in the PROCESSED data ###
##################################

if(knit_manuscript == T){
  render(input = "Manuscript/manuscript_name.Rmd", #Input the path to your .Rmd file here
         output_file = paste("dev_draft_", Sys.Date(), sep = ""), #Name your file here; as it is, this line will create reports named with the date
                                                                  #NOTE: Any file with the dev_ prefix in the Drafts directory will be ignored. Remove "dev_" if you want to include draft files in the GitHub repo
         output_dir = "Output/Drafts/", #Set the path to the desired output directory here
         output_format = "all",
         clean = T)
}
