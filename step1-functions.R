### Author: Selene Banuelos
### Date: 06/01/2023
### Description: Custom functions for cleaning luminex data

# wraps all functions together to reformat data
# timepoint_platform should be a string with the format: 'T2_4plex'
reformatRawData <- function(timepoint_platform, 
                            plasma # dataframe with plasma availability
                            ) {
    #timepoint_platform <- "T4_4plex"
    
    data_list <- importDataAsList(timepoint_platform, "\\.csv") %>%
        lapply(renameCHI3L1) %>%
        lapply(renameIL1B) %>%
        lapply(renameTNFa) %>%
        lapply(changeNaNtoNA) %>%
        lapply(addBatchDateAnalytes) %>%
        lapply(subsetAndCombineData) %>%
        lapply(filter, specimenID != 'blank') %>%
        lapply(addTimePointCol, plasma) %>%
        lapply(rename, dilution_factor = 'Dilution Factor') %>%
        lapply(compareSampleToStd) %>%
        lapply(changeNaNtoNA)
    
    
    #check <- data_list[[1]]
}

# import all data within directory into a list of dataframes - **edit this later to give the option of adding column with file name or not**
importDataAsList <- function(directory, file_extension) {
    # directory <- 'T4_4plex'
    # file_extension <- '\\.csv'
    
    # save .csv file names as list
    file_names <- list.files(path = directory, pattern = file_extension, full.names = TRUE)
    
    # Pick on of the options below (comment out unused line)
    # import all data csv files into list of data frames with column that contains file names
    data_in_list <- lapply(file_names, addColWithFileName)
    
    # # import all data csv files into list of data frames without column that contains file names
    # data_in_list <- lapply(file_names, read.csv)
}

# Function that adds original .csv file name as column in dataframe 
addColWithFileName <- function(x) {
    #x <- file_names[[2]]
    
    df_with_filename <- read.csv(x) %>%
        mutate(file_name = x)
}

# Correct all misspellings of CHI3L1 to 'CHI3L1'
renameCHI3L1 <- function(df) {
    return <- df %>%
        mutate_all(stringr::str_replace_all, "CH.*L1$", "CHI3L1")# matches any string beginning with "CH" and ending with "L1"
}

# function to make capitalization uniform for IL-1B
renameIL1B <- function(df) {
    return <- df %>%
        mutate_all(stringr::str_replace_all, "IL-1b", "IL-1B")# matches any string beginning with "CH" and ending with "L1"
}

# function to make all variations of TNF-a uniform
renameTNFa <- function(df) {
    return <- df %>%
        mutate_all(stringr::str_replace_all, "TNF.*a$", "TNF-a")# matches any string beginning with "CH" and ending with "L1"
}

# change all NaN to NA in dataframe
changeNaNtoNA <- function(df) {
    #df <- data_list[[9]]
    
    df_changed <- df %>% mutate_all(~ifelse(is.nan(.), NA, .))
}

# Function that adds column populated with batch name(plate name) to each data file
# Reference: https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames
addBatchDateAnalytes <- function(df) {
    # find rows that contain batch date and name
    batch_row <- filter(df, Program == "Batch")
    date_row <- filter(df, Program == "Date")
    
    # obtain batch name from xPONENT column from row found above
    batch_name <- batch_row$xPONENT
    date <- date_row$xPONENT
    
    # extract analyte names and put into string
    analyte_names <- extractAnalyteNames(df) %>%
        toString()
    
    # create new df with added columns populated with batch name and date
    batch_added <- df %>% 
        mutate(Batch = batch_name, .before = Program)
    
    date_added <- batch_added %>% 
        mutate(Date = date, .before = Batch)
    
    analytes_added <- date_added %>%
        mutate(analytes = paste0('analytes_', analyte_names))
}

# Function to create vector with analyte names found in RAW output data file
extractAnalyteNames <- function(df) {
    # df <- unmanipulated output file from xPONENT software
    
    # subset rows that contain protein names 
    subset_rows <- subset(df, subset = Program == "Sample")
    
    # convert first row of dataframe to vector
    row_as_vector <- unname(unlist(subset_rows[1,]))
    
    # find index for name of first protein (located after "Sample")
    first_index <- match("Sample", row_as_vector) + 1
    
    # find index for name of last analyte (located before "Total Events", before first NA, or before first empty string, depending on input dataframe format)
    ifelse(!is.na(match("Total Events", row_as_vector)), # if "Total Events" is found in vector
           (last_index <- match("Total Events", row_as_vector) - 1), # then find index of "Total Events"
           ifelse(anyNA(row_as_vector), # if the row contains any NA's
                  (last_index <- first(which(is.na(row_as_vector))) - 1), # then, find index of first NA
                  (last_index <- first(which(row_as_vector == "")) - 1))) # else, find index of first empty string
    
    # store protein names in vector
    names <- row_as_vector[first_index:last_index]
}

# Function that creates single dataframe containing all desired data types
subsetAndCombineData <- function(df) {
    #df <- t2_3plex_list[[7]]
    
    # create a dataframe with specified data type (net MFI) for standards and samples for all analytes
    net_mfi <- df %>%
        changeBlanksToNA() %>%
        subsetDataType("Net MFI") %>%
        janitor::remove_empty("cols") %>%
        janitor::row_to_names(row_number = 2) %>%
        rename(date = 1, batch = 2) %>%
        rename_at( vars( contains("analytes_") ), list( ~paste0("analytes"))) %>%
        mutate_all(~gsub("analytes_", "", .)) %>%
        rename_at( vars( contains(".csv") ), list( ~paste0("file_name"))) %>%
        select(-'Total Events')
    
    
    # create a dataframe with specified data type (expected concentration) of standards for all analytes
    exp_std_conc <- df %>%
        changeBlanksToNA() %>%
        subsetDataType("Standard Expected Concentration") %>%
        janitor::remove_empty("cols") %>%
        janitor::row_to_names(row_number = 2) %>%
        rename(date = 1, batch = 2) %>%
        rename_at( vars( contains("analytes_") ), list( ~paste0("analytes"))) %>%
        mutate_all(~gsub("analytes_", "", .)) %>%
        rename_at( vars( contains(".csv") ), list( ~paste0("file_name"))) %>%
        rename(Sample = "Reagent")
    
    # create a dataframe with specified data type (%CV replicates) of standards and samples for all analytes
    cv <- df %>%
        changeBlanksToNA() %>%
        subsetDataType("%CV Replicates") %>%
        remove_empty("cols") %>%
        row_to_names(row_number = 2) %>%
        rename(date = 1, batch = 2) %>%
        rename_at( vars( contains("analytes_") ), list( ~paste0("analytes"))) %>%
        mutate_all(~gsub("analytes_", "", .)) %>%
        rename_at( vars( contains(".csv") ), list( ~paste0("file_name")))
    
    # create a dataframe with specified data type (dilution factor) of samples for all analytes
    dilution_factor <- df %>%
        changeBlanksToNA() %>%
        subsetDataType("Dilution Factor") %>%
        janitor::remove_empty("cols") %>%
        janitor::row_to_names(row_number = 2) %>%
        rename(date = 1, batch = 2) %>%
        rename_at( vars( contains("analytes_") ), list( ~paste0("analytes"))) %>%
        mutate_all(~gsub("analytes_", "", .)) %>%
        rename_at( vars( contains(".csv") ), list( ~paste0("file_name"))) %>%
        replace(is.na(.), 0) %>%# replace NA in dilution factor column for standards with 1 since they are undiluted
        rename(specimenID = Sample)
    
    # merge all dataframes from above into one, adding data type as suffix for each variable
    merged0 <- lst(net_mfi, exp_std_conc, cv) %>%
        imap(addDfNameAsSuffix) %>%
        reduce(full_join, by = c('date', 'batch', "Sample", 'file_name', 'analytes')) %>% 
        rename(Location = Location_net_mfi, specimenID = Sample) %>%
        full_join(., dilution_factor, by = c('date', 'batch', 'specimenID', 'Location', 'file_name', 'analytes')) %>%
        mutate_at(names(.) %>%
                      .[! . %in% c('date', 'batch', 'specimenID', 'Location', 'file_name', 'analytes')], as.numeric) # convert numbers in df from characters to numeric class (need to do for model fitting). all columns except "Sample" col contain only numbers
    
    # remove background data 
    merged <- merged0 %>%
        filter(specimenID != 'Background0')
    
    return(merged)
    
    #check <- merged0[[1]]
}

# Function that changes all blank cells in dataframe to NA
changeBlanksToNA <- function(df) {
    df %>%
        mutate_all(na_if, "")
}

#function that subsets specified data type from each plate's output file. Input dataframe must have NA in the place of any blank cells
#Reference: https://stackoverflow.com/questions/48275128/how-to-filter-rows-between-two-specific-values
# df = imported raw data
# data_type = data type of interest you want to contain in new data frame
subsetDataType <- function(df, data_type) {
    # get row number of last plate_data$xPONENT=NA before plate_data$xPONENT="Net MFI"
    last_na_before_data_type <- max(which(is.na(df$xPONENT))[which(is.na(df$xPONENT)) < which(df$xPONENT == data_type)])
    
    # get row number of first row that contains NA in columns 1 and 2 after xPONENT=data_type
    first_na_after_data_type <- min(which(is.na(df$Program))[which(is.na(df$Program)) > which(df$xPONENT == data_type)])
    
    # all rows to remove 
    remove_rows <- c(seq(1:last_na_before_data_type), seq(first_na_after_data_type, nrow(df)))
    
    # create new dataframe without rows selected in step above
    new_df <- df[-remove_rows, ]
}

#Reference: https://stackoverflow.com/questions/65152352/suffixes-when-merging-more-than-two-data-frames-with-full-join
addDfNameAsSuffix <- function(x, y) {
    x %>% rename_with(~paste(., y, sep = '_'), -c('Sample', 'date', 'batch', 'file_name', 'analytes'))
}

# Merge visit number to cleaned data by specimen ID. 
# Purpose of this is to produce column that acts as a sanity check that all the data contained in final cleaned data .csv is from the expected time point
addTimePointCol <- function(df, plasma) {
    merged <- plasma %>%
        rename(specimenID = "specimenid", subjectID = "subjectid") %>%
        mutate(specimenID = as.character(specimenID)) %>%
        right_join(., df, by="specimenID")
}

# Filter the number of standard points given analyte and add binary variable that tells whether or not sample net mfi falls above lowest standard
# The final data list is ready for curve fitting/estimating concentration as well as checking for sample reruns
compareSampleToStd <- function(df) {
    #df <- check
    
    # separate analyte names contained in string to vector of names
    analytes <- strsplit(df[1, 'analytes'], ",") %>%
        unlist() %>%
        str_trim() # remove whitespaces
    
    id_vars <- c('subjectID', 'specimenID', 'plasma', 'visitnum', 'date', 'batch', 'Location', 'file_name', 'analytes', 'dilution_factor')
    
    final_data <- sapply(analytes, function(df) df = data.frame(x = df) ) %>%# create list of data frames, each with data for one analyte
        lapply(., subsetAnalyte, df) %>%
        lapply(removeVarsAnalytePrefix) %>%
        lapply(filterStandards) %>%
        lapply(invalidateResultsCV, id_vars) %>%
        lapply(createLowestStdVar) %>%
        lapply(createBelowLowStdBinary) %>%
        lapply(createHighestStdVar) %>%
        lapply(createAboveHighStdBinary) %>%
        lapply(appendVarsAnalytePrefix, id_vars) %>%
        lapply(select, -contains('_analyte')) %>%
        reduce(full_join, by = id_vars)
    
    #check2 <- final_data[[2]]
}

# append given variables with analyte prefix
appendVarsAnalytePrefix <- function(df, not_cols) {
    #df <- analyte_df[[1]]
    
    analyte <- as.character(df[1, 'analyte'])
    
    cols <- names(df)[!names(df) %in% not_cols]
    
    appended_df <- df %>% 
        rename_with(~paste0(analyte, '_', .x), cols)
}

# create new binary variable which:
# equals 1 if sample net mfi is equal to or above net mfi of highest standard in dataframe
# equals 0 if sample net mfi falls below net mfi of highest standard in dataframe
createAboveHighStdBinary <- function(df){
    #df <- final_data[[1]]
    
    # compare sample net mfi to net mfi of highest standard
    df_var_added <- df %>% 
        mutate('above_highest_std' = case_when(net_mfi > highest_std_mfi ~ 1,
                                               net_mfi <= highest_std_mfi ~ 0)) 
}

# create a new variable in dataframe populated with the net mfi value of the highest standard in dataframe
createHighestStdVar <- function(df){
    #df <-  check2
    
    # std_1 <- with(df, grepl('Standard1', specimenID))
    rows_with_stds <- with(df, grepl('Standard', specimenID))
    
    highest_std <- df %>%
        filter(rows_with_stds) %>%
        select(net_mfi) %>%
        max()
    
    df_var_added <- df %>%
        mutate(highest_std_mfi = highest_std)
}

# create new binary variable which:
# equals 1 if sample net mfi is equal to or above net mfi of lowest standard in dataframe
# equals 0 if sample net mfi falls below net mfi of lowest standard in dataframe
createBelowLowStdBinary <- function(df){
    #df <- final_data[[1]]
    
    # compare sample net mfi to net mfi of lowest standard
    df_var_added_1 <- df %>% 
        mutate('below_lowest_std' = case_when(net_mfi >= lowest_std_mfi ~ 0,
                                              net_mfi < lowest_std_mfi ~ 1,
                                              is.na(lowest_std_mfi) ~ NA))
}

# create a new variable in dataframe populated with the net mfi value of the lowest standard in dataframe
createLowestStdVar <- function(df){
    #df <-  final_data[[1]]
    
    rows_with_stds <- with(df, grepl('Standard', specimenID))
    
    lowest_std <- df %>%
        changeNaNtoNA() %>%
        filter(rows_with_stds) %>%
        select('net_mfi') %>%
        tidyr::drop_na() %>%
        min()
    
    df_var_added <- df %>%
        # mutate(lowest_std_mfi = ifelse(lowest_std >= 10, lowest_std, 'lowest std falls below 10 net mfi')) %>%
        mutate(notes = ifelse(lowest_std < 10, 
                              paste0(df[1, 'analyte'], " lowest standard falls below 10 net mfi"), 
                              NA)) %>%
        mutate(lowest_std_mfi = ifelse(lowest_std >= 10,
                                       lowest_std, NA))
}

# if cv is NA or cv > 20 then invalidate results - one analyte dataframe should be passed in. Also need to specify identification variables so they don't get invalidated.
invalidateResultsCV <- function(one_analyte_df, id_vars) {
    result_vars <- setdiff(names(one_analyte_df), c(id_vars, 'cv', 'analyte'))
    
    check <- one_analyte_df %>%
        mutate(across(result_vars, ~ case_when(cv > 20 ~ NA, 
                                               is.na(cv) ~ NA,
                                               TRUE ~ .)))
}

# THIS FUNCTION IS GROSS, FIX THIS WHEN YOU HAVE TIME
# check out: https://stackoverflow.com/questions/55059761/r-mutate-multiple-columns-with-ifelse-condition
# Function that select correct number of standard points to include in curve, given analyte name
filterStandards <- function(df) {
    stds_5_6_7 <- c("Standard5", "Standard6", "Standard7")# keep standards 1-4
    stds_6_7 <- c("Standard6", "Standard7")# keep standards 1-5
    std_7 <- c("Standard7")# keep standards 1-6
    
    analyte <- df[1, 'analyte']
    
    filtered_stds <- if (analyte == 'CHI3L1'){
        filter(df, (!specimenID %in% stds_6_7))
    } 
    else if (analyte == 'BDNF') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'MPO') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'CRP') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'IGFBP-1') {
        filter(df, (!specimenID %in% stds_6_7))
    }
    else if (analyte == 'IGFBP-3') {
        filter(df, (!specimenID %in% stds_5_6_7))
    }
    else if (analyte == 'Leptin') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'TNF-a') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'VEGF') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'IL-1B') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'IL-6') {
        filter(df, (!specimenID %in% std_7))
    }
    else if (analyte == 'IL-10') {
        filter(df, (!specimenID %in% std_7))
    }
    else {
        df
    }
}

# remove analyte prefix from vars to make manipulations easier
removeVarsAnalytePrefix <- function(df) {
    #df <- analyte_df[[1]]
    analyte <- as.character(df[1, 'analyte'])
    
    prefix <- paste0(analyte, '_')
    
    remove_df <- df %>%
        rename_all(~str_replace(., paste0(analyte, '_'), ""))
}

# function that creates separate data frames for each analyte
subsetAnalyte <- function(df, plate_data) {
    analyte_name <- df
    
    analyte_df <- plate_data %>%
        select(c('subjectID', 'specimenID', 'plasma', 'visitnum', 'date', 'batch', 'Location', 'dilution_factor', 'file_name', 'analytes', contains(analyte_name))) %>%
        #rename_at(vars(contains(analyte_name)), ~ str_replace(., paste0(analyte_name, "_"), "")) %>%
        mutate(analyte = analyte_name)
}