### Author: Selene Banuelos
### Date: 06/01/2023
### Description: step 1 of Luminex data cleaning pipeline

# step 1 -----------------------------------------------------------------------
# setup
library(dplyr)
library(purrr)
library(janitor)
library(stringr)
library(tidyr)
source('code/step1-functions.R') # import helper functions

# function used to reformat raw xPONENT data
# path should be a string with the format: 'T2_4plex'
reformatRawData <- function(path, 
                            plasma # dataframe with plasma availability
                            ){

    # import all csv in path as data frames into a list
    data_list <- importDataAsList(path, "\\.csv") %>%
        # standardize spelling of analytes
        lapply(function(df) df %>% 
                   # correct all misspellings of CHI3L1 to 'CHI3L1'
                   mutate_all(str_replace_all, "CH.*L1$", "CHI3L1") %>%
                   # make capitalization uniform for IL-1B
                   mutate_all(str_replace_all, "IL-1b", "IL-1B") %>%
                   # make all variations of TNF-a uniform
                   mutate_all(str_replace_all, "TNF.*a$", "TNF-a")
        ) %>%
        # change all NaN to NA in dataframe
        lapply(function(df) mutate_all(df, ~ifelse(is.nan(.), NA, .))) %>%
        # extract batch name, run date, and analyte information
        lapply(extractBatchDateAnalytes) %>%
        # extract net mfi, %CV, std concentration, and dilution factor data
        lapply(extractAndCombineData) %>% 
        # remove empty wells
        lapply(filter, specimenID != 'blank') %>%
        # Merge visit number to cleaned data (do i need this here?)
        lapply(addVisitnum, plasma) %>%
        # do i need to do this?
        lapply(rename, dilution_factor = 'Dilution Factor')
    
}

# import data
# plasma availability data
plasma <- read.csv('data-raw/PEARLSBio-Plasma_DATA_2023-06-28_1451.csv')

# extract necessary data from raw files & format for standard curve fitting and
# estimating concentrations

# timepoint 2
t2_3plex_list <- reformatRawData('data-raw/T2_3plex', plasma)
t2_4plex_list <- reformatRawData('data-raw/T2_4plex', plasma)
t2_6plex_list <- reformatRawData('data-raw/T2_6plex', plasma)

# timepoint 4
t4_3plex_list <- reformatRawData('data-raw/T4_3plex', plasma)
t4_4plex_list <- reformatRawData('data-raw/T4_4plex', plasma)
t4_hscytokine_list <- reformatRawData('data-raw/T4_HScytokine', plasma)

# remove unwanted objects from environment
rm(list = setdiff(ls(), ls(pattern = '_list')))

# save data objects for next step
save.image('data-processed/step1-data.RData')