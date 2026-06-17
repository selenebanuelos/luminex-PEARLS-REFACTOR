### Author: Selene Banuelos
### Date: 06/01/2023
### Description: Custom functions for cleaning luminex data

# wraps all functions together to reformat data
# timepoint_platform should be a string with the format: 'T2_4plex'
reformatRawData <- function(path, 
                            plasma # dataframe with plasma availability
                            ) {
    # # for testing
    # path <- "data-raw/T2_6plex"
    # plasma <- plasma
    
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
    
    #check <- data_list[[1]]
}

# import all data within directory into a list of data frames with file names
importDataAsList <- function(directory, file_extension) {
    
    # for testing
    # directory <- 'data-raw/T2_3plex'
    # file_extension <- '\\.csv'
    
    # save .csv file names as list
    file_names <- list.files(path = directory, pattern = file_extension, full.names = TRUE)
    
    # import all data csv files into list of data frames
    data_in_list <- lapply(file_names,
                           # add original .csv file name as column in dataframe 
                           function(x) read.csv(x) %>% mutate(file_name = x)
                           )
}

# Function that adds column populated with batch name(plate name) to each data file
extractBatchDateAnalytes <- function(df) {
    
    # # for testing
    # df <- check
    
    # extract batch name and date of run as strings (values in 'xPONENT' column)
    batch_name <- filter(df, Program == "Batch") %>% pull(xPONENT)
    date_run <- filter(df, Program == "Date") %>% pull(xPONENT)
    
    # create string that contains analyte names, separated by commas
    analyte_names <- df %>%
        # analyte names stored in rows where column 'Program' == 'Sample'
        filter(Program == 'Sample') %>%
        # remove row name ('Sample') and file name
        select(-c(Program, file_name)) %>%
        # change any empty strings to NA
        mutate_all(list(~ na_if(., ''))) %>%
        # remove empty NA columns
        remove_empty('cols') %>%
        # keep first row
        slice(1) %>%
        # convert row into vector
        unlist(.) %>%
        # combine analyte names in vector into one string, separated by commas
        toString(.)
    
    # return df with added columns describing batch, date, and analyte names
    df %>% 
        # add in descriptive suffix for downstream data manipulation
        mutate(batch = batch_name,
               date = date_run,
               analytes = analyte_names
               )
}

# Function that creates single dataframe containing all desired data types
extractAndCombineData <- function(df) {
    
    # # for testing
    # df <- check
    
    # column names of file name, batch, date, analyte data
    info_vars <- c('file_name', 'batch', 'date', 'analytes')
    
    # separate run info from data
    info <- select(df, all_of(info_vars))
    data <- select(df, -info_vars) %>% 
        # change any empty strings to NA
        mutate_all(na_if, '') 
    
    # function that subsets specified data type from each plate's raw file 
    # Input dataframe must have NA in the place of any blank cells
    extractDataType <- function(df, # data frame
                                data_type # data type as string
                                ) {
        # get row number of last plate_data$xPONENT=NA before plate_data$xPONENT="Net MFI"
        last_na_before_data_type <- max(which(is.na(df$xPONENT))[which(is.na(df$xPONENT)) < which(df$xPONENT == data_type)])
        
        # get row number of first row that contains NA in columns 1 and 2 after xPONENT=data_type
        first_na_after_data_type <- min(which(is.na(df$Program))[which(is.na(df$Program)) > which(df$xPONENT == data_type)])
        
        # all rows to remove 
        remove_rows <- c(seq(1:last_na_before_data_type), 
                         seq(first_na_after_data_type, nrow(df))
                         )
        
        # new dataframe without rows selected in step above
        df[-remove_rows, ] %>%
            # remove columns with no data
            remove_empty("cols") %>%
            # use second row as variable names (first row just specifies data type)
            row_to_names(row_number = 2)
            
    }
    
    # create a dataframe with specified data type (net MFI) for standards and samples for all analytes
    net_mfi <- data %>%
        # extract rows that contain Net MFI-related data
        extractDataType("Net MFI") %>%
        # remove unecessary Total Events column
        select(-'Total Events')
    
    # create a dataframe with specified data type (expected concentration) of standards for all analytes
    exp_std_conc <- data %>%
        extractDataType('Standard Expected Concentration') %>%
        # remove unecessary Reagent column
        rename(Sample = "Reagent")
    
    # create a dataframe with specified data type (%CV replicates) of standards and samples for all analytes
    cv <- data %>%
        extractDataType('%CV Replicates')
    
    # create a dataframe with specified data type (dilution factor) of samples for all analytes
    dilution_factor <- data %>%
        extractDataType("Dilution Factor") %>%
        # replace NA in dilution factor column for standards with 1 since they are undiluted
        replace(is.na(.), 0) %>% 
        # rename specimenID for donwstream joining
        rename(specimenID = Sample)
    
    # merge all dataframes from above into one
    merged <- lst(net_mfi, exp_std_conc, cv) %>%
        # add data frame names as suffix to all variables except 'Sample'
        imap(function(df, df_name) rename_with(df, 
                                               ~paste(., df_name, sep = '_'), 
                                               -c('Sample'))) %>%
        # join net mfi, standards, and %CV data together
        reduce(full_join, by = 'Sample') %>% 
        # rename location and sample variables for downstream joining
        rename(Location = Location_net_mfi, 
               specimenID = Sample) %>%
        # join with dilution factor data
        full_join(dilution_factor, 
                  by = c('specimenID', 'Location')) %>%
        # change specific variables into numeric class
        mutate_at(names(.) %>%
                      .[! . %in% c('specimenID', 'Location')], 
                  as.numeric) %>%
        # remove background data
        filter(specimenID != 'Background0') %>%
        # add in run info back into data
        cbind(., slice(info, 1:nrow(.)))
    
    return(merged)
    
    #check <- merged0[[1]]
}

# Purpose of this is to add column that acts as a sanity check that all the 
# data contained in final cleaned data .csv is from the expected time point
addVisitnum <- function(df, plasma) {
    merged <- plasma %>%
        rename(specimenID = "specimenid", subjectID = "subjectid") %>%
        mutate(specimenID = as.character(specimenID)) %>%
        right_join(., df, by="specimenID")
}

# # Filter the number of standard points given analyte and add binary variable that tells whether or not sample net mfi falls above lowest standard
# # The final data list is ready for curve fitting/estimating concentration as well as checking for sample reruns
# compareSampleToStd <- function(df) {
#     #df <- check
#     
#     # separate analyte names contained in string to vector of names
#     analytes <- strsplit(df[1, 'analytes'], ",") %>%
#         unlist() %>%
#         str_trim() # remove whitespaces
#     
#     id_vars <- c('subjectID', 'specimenID', 'plasma', 'visitnum', 'date', 'batch', 'Location', 'file_name', 'analytes', 'dilution_factor')
#     
#     final_data <- sapply(analytes, function(df) df = data.frame(x = df) ) %>%# create list of data frames, each with data for one analyte
#         lapply(., subsetAnalyte, df) %>%
#         lapply(removeVarsAnalytePrefix) %>%
#         lapply(filterStandards) %>%
#         lapply(invalidateResultsCV, id_vars) %>%
#         lapply(createLowestStdVar) %>%
#         lapply(createBelowLowStdBinary) %>%
#         lapply(createHighestStdVar) %>%
#         lapply(createAboveHighStdBinary) %>%
#         lapply(appendVarsAnalytePrefix, id_vars) %>%
#         lapply(select, -contains('_analyte')) %>%
#         reduce(full_join, by = id_vars)
#     
#     #check2 <- final_data[[2]]
# }
# 
# # append given variables with analyte prefix
# appendVarsAnalytePrefix <- function(df, not_cols) {
#     #df <- analyte_df[[1]]
#     
#     analyte <- as.character(df[1, 'analyte'])
#     
#     cols <- names(df)[!names(df) %in% not_cols]
#     
#     appended_df <- df %>% 
#         rename_with(~paste0(analyte, '_', .x), cols)
# }
# 
# # create new binary variable which:
# # equals 1 if sample net mfi is equal to or above net mfi of highest standard in dataframe
# # equals 0 if sample net mfi falls below net mfi of highest standard in dataframe
# createAboveHighStdBinary <- function(df){
#     #df <- final_data[[1]]
#     
#     # compare sample net mfi to net mfi of highest standard
#     df_var_added <- df %>% 
#         mutate('above_highest_std' = case_when(net_mfi > highest_std_mfi ~ 1,
#                                                net_mfi <= highest_std_mfi ~ 0)) 
# }
# 
# # create a new variable in dataframe populated with the net mfi value of the highest standard in dataframe
# createHighestStdVar <- function(df){
#     #df <-  check2
#     
#     # std_1 <- with(df, grepl('Standard1', specimenID))
#     rows_with_stds <- with(df, grepl('Standard', specimenID))
#     
#     highest_std <- df %>%
#         filter(rows_with_stds) %>%
#         select(net_mfi) %>%
#         max()
#     
#     df_var_added <- df %>%
#         mutate(highest_std_mfi = highest_std)
# }
# 
# # create new binary variable which:
# # equals 1 if sample net mfi is equal to or above net mfi of lowest standard in dataframe
# # equals 0 if sample net mfi falls below net mfi of lowest standard in dataframe
# createBelowLowStdBinary <- function(df){
#     #df <- final_data[[1]]
#     
#     # compare sample net mfi to net mfi of lowest standard
#     df_var_added_1 <- df %>% 
#         mutate('below_lowest_std' = case_when(net_mfi >= lowest_std_mfi ~ 0,
#                                               net_mfi < lowest_std_mfi ~ 1,
#                                               is.na(lowest_std_mfi) ~ NA))
# }
# 
# # create a new variable in dataframe populated with the net mfi value of the lowest standard in dataframe
# createLowestStdVar <- function(df){
#     #df <-  final_data[[1]]
#     
#     rows_with_stds <- with(df, grepl('Standard', specimenID))
#     
#     lowest_std <- df %>%
#         mutate_all(~ifelse(is.nan(.), NA, .)) %>%
#         filter(rows_with_stds) %>%
#         select('net_mfi') %>%
#         tidyr::drop_na() %>%
#         min()
#     
#     df_var_added <- df %>%
#         # mutate(lowest_std_mfi = ifelse(lowest_std >= 10, lowest_std, 'lowest std falls below 10 net mfi')) %>%
#         mutate(notes = ifelse(lowest_std < 10, 
#                               paste0(df[1, 'analyte'], " lowest standard falls below 10 net mfi"), 
#                               NA)) %>%
#         mutate(lowest_std_mfi = ifelse(lowest_std >= 10,
#                                        lowest_std, NA))
# }
# 
# # if cv is NA or cv > 20 then invalidate results - one analyte dataframe should be passed in. Also need to specify identification variables so they don't get invalidated.
# invalidateResultsCV <- function(one_analyte_df, id_vars) {
#     result_vars <- setdiff(names(one_analyte_df), c(id_vars, 'cv', 'analyte'))
#     
#     check <- one_analyte_df %>%
#         mutate(across(result_vars, ~ case_when(cv > 20 ~ NA, 
#                                                is.na(cv) ~ NA,
#                                                TRUE ~ .)))
# }
# 
# # THIS FUNCTION IS GROSS, FIX THIS WHEN YOU HAVE TIME
# # check out: https://stackoverflow.com/questions/55059761/r-mutate-multiple-columns-with-ifelse-condition
# # Function that select correct number of standard points to include in curve, given analyte name
# filterStandards <- function(df) {
#     stds_5_6_7 <- c("Standard5", "Standard6", "Standard7")# keep standards 1-4
#     stds_6_7 <- c("Standard6", "Standard7")# keep standards 1-5
#     std_7 <- c("Standard7")# keep standards 1-6
#     
#     analyte <- df[1, 'analyte']
#     
#     filtered_stds <- if (analyte == 'CHI3L1'){
#         filter(df, (!specimenID %in% stds_6_7))
#     } 
#     else if (analyte == 'BDNF') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'MPO') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'CRP') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'IGFBP-1') {
#         filter(df, (!specimenID %in% stds_6_7))
#     }
#     else if (analyte == 'IGFBP-3') {
#         filter(df, (!specimenID %in% stds_5_6_7))
#     }
#     else if (analyte == 'Leptin') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'TNF-a') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'VEGF') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'IL-1B') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'IL-6') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else if (analyte == 'IL-10') {
#         filter(df, (!specimenID %in% std_7))
#     }
#     else {
#         df
#     }
# }
# 
# # remove analyte prefix from vars to make manipulations easier
# removeVarsAnalytePrefix <- function(df) {
#     #df <- analyte_df[[1]]
#     analyte <- as.character(df[1, 'analyte'])
#     
#     prefix <- paste0(analyte, '_')
#     
#     remove_df <- df %>%
#         rename_all(~str_replace(., paste0(analyte, '_'), ""))
# }
# 
# # function that creates separate data frames for each analyte
# subsetAnalyte <- function(df, plate_data) {
#     analyte_name <- df
#     
#     analyte_df <- plate_data %>%
#         select(c('subjectID', 'specimenID', 'plasma', 'visitnum', 'date', 'batch', 'Location', 'dilution_factor', 'file_name', 'analytes', contains(analyte_name))) %>%
#         #rename_at(vars(contains(analyte_name)), ~ str_replace(., paste0(analyte_name, "_"), "")) %>%
#         mutate(analyte = analyte_name)
# }

### REMOVE THESE FUNCTIONS #####################################################
# Function that adds original .csv file name as column in dataframe 
# addColWithFileName <- function(x) {
#     #x <- file_names[[2]]
#     
#     df_with_filename <- read.csv(x) %>%
#         mutate(file_name = x)
# }

# # Correct all misspellings of CHI3L1 to 'CHI3L1'
# renameCHI3L1 <- function(df) {
#     return <- df %>%
#         mutate_all(stringr::str_replace_all, "CH.*L1$", "CHI3L1")# matches any string beginning with "CH" and ending with "L1"
# }

# # function to make capitalization uniform for IL-1B
# renameIL1B <- function(df) {
#     return <- df %>%
#         mutate_all(stringr::str_replace_all, "IL-1b", "IL-1B")# matches any string beginning with "CH" and ending with "L1"
# }

# # function to make all variations of TNF-a uniform
# renameTNFa <- function(df) {
#     return <- df %>%
#         mutate_all(stringr::str_replace_all, "TNF.*a$", "TNF-a")# matches any string beginning with "CH" and ending with "L1"
# }

# # change all NaN to NA in dataframe
# changeNaNtoNA <- function(df) {
#     #df <- data_list[[9]]
#     
#     df_changed <- df %>% mutate_all(~ifelse(is.nan(.), NA, .))
# }

# # Function that changes all blank cells in dataframe to NA
# changeBlanksToNA <- function(df) {
#     df %>%
#         mutate_all(na_if, "")
# }

# # Function to create vector with analyte names found in RAW output data file
# extractAnalyteNames <- function(df) {
#     
#     # for testing
#     #df <- check # unmanipulated output file from xPONENT software
#     
#     # create string that contains analyte names, separated by commas
#     df %>%
#         # analyte names stored in rows where column 'Program' == 'Sample'
#         filter(Program == 'Sample') %>%
#         # remove row name ('Sample') and file name
#         select(-c(Program, file_name)) %>%
#         # change any empty strings to NA
#         mutate_all(list(~ na_if(., ''))) %>%
#         # remove empty NA columns
#         remove_empty('cols') %>%
#         # keep first row
#         slice(1) %>%
#         # convert row into vector
#         unlist(.) %>%
#         # combine analyte names in vector into one string, separated by commas
#         toString(.)
# }
