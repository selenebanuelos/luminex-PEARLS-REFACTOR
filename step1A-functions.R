### Author: Selene Banuelos
### Date: 06/01/2023
### Description: Custom functions for cleaning luminex data

# NEED TO ADD IN: 
#     
#     lapply(compareSampleToStd) %>%
#     # change all NaN to NA in dataframe
#     lapply(function(df) mutate_all(df, ~ifelse(is.nan(.), NA, .)))

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
        mutate_all(~ifelse(is.nan(.), NA, .)) %>%
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