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