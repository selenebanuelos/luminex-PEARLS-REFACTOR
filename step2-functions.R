### Author: Selene Banuelos
### Date: 06/01/2023
### Description: Custom functions for cleaning luminex data

# function that warps all functions together: 
# takes in a single plate's data and fits a standard curve to produce estimates 
# of concentrations
estimateConcentrations <- function(df){
    #df <- t2_6plex_list[[1]]
    
    # separate analyte names contained in string to vector of names
    analytes <- strsplit(df$analytes[[1]], ",") %>%
        unlist() %>%
        str_trim() # remove whitespaces
    
    id_vars <- c('date', 'batch', 'Location', 'plasma', 'subjectID', 'specimenID', 'visitnum', 'dilution_factor', 'file_name', 'analytes')
    
    # create a dataframe with net MFI and %CV for standards and samples for all analytes as well as expected concentration for all standards
    final_data <- sapply(analytes, function(df_0) df_0 = data.frame(x = df_0) ) %>%# create list of data frames, each with data for one analyte
        lapply(subsetAnalyte, df) %>%
        lapply(removeVarsAnalytePrefix) %>%
        lapply(filterStandards) %>%
        lapply(fitCurveEstimateConc) %>%
        lapply(createLowestStdConc) %>%
        lapply(createHighestStdConc) %>%
        lapply(appendVarsAnalytePrefix, id_vars) %>%
        lapply(select, -contains('_analyte')) %>%
        reduce(full_join, by = id_vars)
}

# function that creates separate data frames for each analyte
subsetAnalyte <- function(df, plate_data) {
    analyte_name <- df
    
    analyte_df <- plate_data %>%
        select(c('subjectID', 'specimenID', 'plasma', 'visitnum', 'date', 'batch', 'Location', 'dilution_factor', 'file_name', 'analytes', contains(analyte_name))) %>%
        #rename_at(vars(contains(analyte_name)), ~ str_replace(., paste0(analyte_name, "_"), "")) %>%
        mutate(analyte = analyte_name)
}

# remove analyte prefix from vars to make manipulations easier
removeVarsAnalytePrefix <- function(df) {
    #df <- analyte_df[[1]]
    analyte <- as.character(df[1, 'analyte'])
    
    prefix <- paste0(analyte, '_')
    
    remove_df <- df %>%
        rename_all(~str_replace(., paste0(analyte, '_'), ""))
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

# fit standard curve and estimate concentrations of unknown specimen + standard point concentrations
# cannot input NA values into getEstimates, will generate error
fitCurveEstimateConc <- function(df){
    #df <- final_data[[1]]
    
    analyte <- df[1, 'analyte']
    
    # keep specimen with net_mfi NA in separate data frame
    df_NA <- df %>% filter(is.na(net_mfi))
    
    df_for_model <- df %>% 
        filter(!is.na(net_mfi)) %>%
        generateDfForModel()
    
    # generate model to fit standard curve based on standard points
    model <- generateModel(df_for_model)
    
    # nplr::getEstimates() inverts the function and returns the estimated concentration for a given response (MFI value)
    # estimate concentrations of samples/standards using mfi values and standard curve
    # cannot input NA values into getEstimates, will generate error
    estimated_conc <- getEstimates(model, df_for_model$net_mfi_prop, conf.level=0.95) %>%
        rename(net_mfi_prop_used = "y", estimated_concentration = "x")
    
    # combine estimated concentrations to original data frame (only keep rows with "net_mfi_prop" values found in both df)
    merged <- bind_cols(df_for_model, estimated_conc) %>% # net_mfi_prop should match net_mfi_prop_used. These were included as sanity checks
        mutate_at(names(.) %>%
                      .[! . %in% c('date', 'batch', 'Location','subjectID', 'specimenID', paste0(analyte, '_notes'), 'file_name', 'analytes', 'analyte')], as.numeric) %>%# convert numbers in df from characters to numeric class (need to do for model fitting). all columns except "Sample" col contain only numbers
        mutate(est_conc_times_dilution_factor = estimated_concentration * (dilution_factor + 1)) %>% # dilution factor of 1 means 1:1 dilution
        rename(CI_lower_bound = x.025, CI_upper_bound = x.975)
    
    final_data_0 <- full_join(merged, df_NA, by = names(df_NA))
    
}

# FUNCTIONS FROM CURVE FITTING AND ESTIMATE CONCENTRATIONS SCRIPT ##############

# format df for nplr(): repeating x variable (exp standard conc) in first column, resulting MFI in second column (y)
generateDfForModel <- function(df){
    #df <-  df_for_model
    
    # identify rows with data pertaining to standards
    rows_with_stds <- with(df, grepl("Standard", specimenID))
    
    # identify max mfi value of standards to use as max value when converting mfi to proportions
    max_mfi <- df %>%
        filter(rows_with_stds) %>%
        select(contains("net_mfi")) %>%
        max()
    
    # create variable with all mfi values converted to proportions using nplr::convertToProp()
    df_for_model <- df %>%
        mutate(net_mfi_prop = convertToProp(select(., contains("net_mfi")), T0 = 0, Ctrl = max_mfi)) %>%
        mutate(net_mfi_prop = pull(pull(., "net_mfi_prop"))) # need to convert this variable to vector since nplr::convertToProp() created it as data frame within data frame
}

# fit standard curve with nplr() 
generateModel <- function(df){
    #df <- df_for_model
    
    # identify rows with data pertaining to standards
    rows_with_stds <- with(df, grepl("Standard", specimenID))
    
    # x = expected standard concentrations
    x_exp_std_conc_0 <- df %>%
        select(exp_std_conc) %>%
        .[rows_with_stds, ]# keep only rows pertaining to standards and columns with exp_conc and net_mfi data
    
    # check that 'x_exp_std_conc_0' is vector, if not, turn it into a vector using pull()
    # some weird, inconsistent behavior was happening so I needed to add this step. This is a bandaid and it would be ideal to figure out why a df is created some times and a vector other times
    if(is.data.frame(x_exp_std_conc_0)){
        x_exp_std_conc <- pull(x_exp_std_conc_0, exp_std_conc)
    } else {
        x_exp_std_conc <- x_exp_std_conc_0
    }
    
    # y = average net MFI values for each standard. MFI values must be converted to proportions.
    y_std_net_mfi_0 <- df %>%
        select(net_mfi_prop) %>%
        .[rows_with_stds, ] # keep only rows pertaining to standards and columns with exp_conc and net_mfi data
    
    # check that 'y_std_net_mfi_0' is vector, if not, turn it into a vector using pull()
    # some weird, inconsistent behavior was happening so I needed to add this step. This is a bandaid and it would be ideal to figure out why a df is created some times and a vector other times
    if(is.data.frame(y_std_net_mfi_0)){
        y_std_net_mfi <- pull(y_std_net_mfi_0, net_mfi_prop)
    } else {
        y_std_net_mfi <- y_std_net_mfi_0
    }
    
    # print analyte name
    print('#####################################################################')
    print(paste0(df[1, 'file_name'], ' ', df[1, 'analyte'], ' Standard Curve:'))
    
    # find best model using nplr()
    model <- nplr(x_exp_std_conc, y_std_net_mfi)
    
    # have to call the model to fit curve
    print(model)
}

# create a new variable in dataframe populated with the estimated concentration of the lowest standard in dataframe
createLowestStdConc <- function(df){
    #df <-  final_data[[1]]
    
    rows_with_stds <- with(df, grepl('Standard', specimenID))
    
    lowest_std <- df %>%
        changeNaNtoNA() %>%
        filter(rows_with_stds) %>%
        select(estimated_concentration) %>%
        drop_na() %>%
        min()
    
    df_var_added <- df %>%
        mutate(lowest_std_conc = lowest_std)
}

# create a new variable in dataframe populated with the estimated concentration of the highest standard in dataframe
createHighestStdConc <- function(df){
    #df <-  final_data[[1]]
    
    rows_with_stds <- with(df, grepl('Standard', specimenID))
    
    highest_std <- df %>%
        filter(rows_with_stds) %>%
        select(estimated_concentration) %>%
        max()
    
    df_var_added <- df %>%
        mutate(highest_std_conc = highest_std)
}

# append given variables with analyte prefix
appendVarsAnalytePrefix <- function(df, not_cols) {
    #df <- analyte_df[[1]]
    
    analyte <- as.character(df[1, 'analyte'])
    
    cols <- names(df)[!names(df) %in% not_cols]
    
    appended_df <- df %>% 
        rename_with(~paste0(analyte, '_', .x), cols)
}