### Author: Selene Banuelos
### Date: 06/01/2023
### Description: 

# step 1
################################################################################
# setup
library(dplyr)
library(purrr)
library(janitor)
library(stringr)
library(tidyr)
source('code/step1-functions.R')

# import data
# plasma availability data
plasma <- read.csv('data-raw/PEARLSBio-Plasma_DATA_2023-06-28_1451.csv')

# extract necessary data from raw files & format for standard curve fitting and
# estimating concentrations
t2_3plex_list <- reformatRawData('data-raw/T2_3plex', plasma)

# remove unwanted objects from environment
rm(list = setdiff(ls(), ls(pattern = '_list')))

# save data objects for next step
save.image('data-processed/step1-data.RData')

# step 2
################################################################################
# setup
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(nplr)
source('code/step2-functions.R')

# load data from previous step
load('data-processed/step1-data.RData')

# Using full data (with standards and controls), fit standard curve and estimate 
# concentrations of both unknown specimen and standard points. Recombine dataframes
# (according to timepoint and assay) from manually estimated concentrations 
# ('problem' data) and automated estimated concentrations
t2_3plex_estimated <- map_df(t2_3plex_list, estimateConcentrations)

# remove unwanted objects from environment
rm(list = setdiff(ls(), ls(pattern = '_estimated')))

# save data objects for next step
save.image('data-processed/step2-data.RData')

# step X - dont think this is necessary
################################################################################
# Double check that dilution factors look reasonable before selecting reruns
# define function
checkDilFactor <- function(df) {
    check <- df %>%
        select('dilution_factor') %>%
        unique()
}

t2_3_df <- checkDilFactor(t2_3plex_estimated)

# step 3
################################################################################
# setup
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
source('code/step3-functions.R')

# load data from previous step
load('data-processed/step2-data.RData')

# For samples that were run in more than 1 plate, identify which rerun to keep 
# based on algorithm shown in flowchart
# timepoint 2
t2_3plex_reruns_annotated <- selectRerunToKeep(t2_3plex_estimated)

# remove unwanted objects from environment
rm(list = setdiff(ls(), ls(pattern = '_reruns_annotated')))

# save data objects for next step
save.image('data-processed/step3-data.RData')

# step 4
################################################################################
# setup
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
source('code/step4-functions.R')

# load data from previous step
load('data-processed/step3-data.RData')

# At this step, you may be mixing analytes from different batches for one 
# specimen ID. For example if run on date 1/2/2020 had QC-passing data for CRP
# and MPO for specimenID 500 (but not BDNF), the BDNF data (if passess QC) in 
# specimen 500 from date 3/4/2020 will be appended to the earlier date. Date, 
# Location, batch, file_name, dilution factor variables will be dropped at this 
# point to allow for combinations of data from different run dates.
t2_3plex_cleaned <- cleanData(t2_3plex_reruns_annotated)

# remove unwanted objects from environment
rm(list = setdiff(ls(), ls(pattern = '_cleaned')))

# save data objects for next step
save.image('data-processed/step4-data.RData')