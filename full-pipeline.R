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
source('code/step1-functions.R')

# import data
# plasma availability data
plasma <- read.csv('data-raw/PEARLSBio-Plasma_DATA_2023-06-28_1451.csv')

# extract necessary data from raw files & format for standard curve fitting and
# estimating concentrations
t2_3plex_list <- reformatRawData('data-raw/T2_3plex', plasma)

# save data objects for next step
save.image('data-processed/step1-data.RData')

# step 2
################################################################################
# setup
library(dplyr)
library(nplr)
source('code/step2-functions.R')

# load data from previous step
load('data-processed/step1-data.RData')

# Using full data (with standards and controls), fit standard curve and estimate 
# concentrations of both unknown specimen and standard points. Recombine dataframes
# (according to timepoint and assay) from manually estimated concentrations 
# ('problem' data) and automated estimated concentrations
t2_3plex_estimated <- map_df(t2_3plex_list, estimateConcentrations)

# step 3
################################################################################
# Double check that dilution factors look reasonable before selecting reruns
# define function
checkDilFactor <- function(df) {
    check <- df %>%
        select('dilution_factor') %>%
        unique()
}

t2_3_df <- checkDilFactor(t2_3plex_estimated)