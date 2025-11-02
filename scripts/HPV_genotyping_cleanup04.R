# Preamble ---------------------------------------------------------------------

# Code author   : Khomotso Mohlala
# Project       : Patterns of HPV co-infections, persistence, and infection 
#                 trajectories in South African women
# Supervisors   : Prof. Lenine Liebenberg, Dr. Eduan Wilkinson, 
#                 Dr. Tomasz Sanko, Dr. Cari van Schalkwyk
# Description   : This script contains reproducible code for cleaning data 
#                 collected as part of the FRESH cohort in KwaZulu-Natal (KZN).


# Notes:
#  - Only records collected before 2022 are retained, as the project investigates
#    patterns of HPV before the vaccine roll-out. Participants enrolled after
#    this point likely received HPV vaccines.
#  - Data sets are imported from structured Excel files.
#  - This script handles the initial data cleaning processes.

# Setup ------------------------------------------------------------------------
# Libraries used. If not installed, use: install.packages("package_name")

library(readxl)      #Read Excel files
library(tidyverse)   #Tidy data tools
library(janitor)     #Clean data frames
library(readr)       #Read and write CSV files

#-------------------------------------------------------------------------------
#set working directory to folder that contains excel files

WKDIR <- "C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts"
setwd(WKDIR)

#-------------------------------------------------------------------------------
#load data into R using excel file name and view table  

filename <- "04. HPV_genotyping_data.xlsx"

# Check if file exists before attempting to read it

if (file.exists(filename)) {
  data_04  <- read_excel(filename, sheet = 'curated prevax_2019specimens')
  View(data_04)
} else {
  stop("File not found: ", filename)
}

#------------------------------------------------------------------------------- 
# Keep only key variables for cleaning

HPV_genotyping_data <- data_04 %>% 
  select(
      study_id,
      visit_code,
      large_blood_date,
      nugent_gardnerella,
      nugent_mobiluncus,
      nugent_lactobacillus,
      nugent_total_score,
      nugent_interpretation,
      sti_mycoplasma_genitalium_pcr,
      sti_trichomonas_vaginalis_pcr,
      sti_neisseria_gonorrhoeae_pcr,
      sti_chlamydia_trachomatis_pcr,
      # NOTE: Excluded HPV types 44,55 62 and 81
      # Genotyping assay used cannot discriminate between pairs therefore 
      # Genotypes may represent either single or multiple infections due to this uncertainty
      # we exclude them from ALL analysis 
      `16`, `18`, `26`, `31`, `33`, `35`, `39`, `45`, `51`, `52`, `53`,
      `56`, `58`, `59`, `66`, `68`, `73`, `82`, `6`, `11`, `40`, `42`,
      `43`, `54`, `61`, `67`, `69`, `70`, `71`, `72`, `84`,
      REPORT,
      COMMENT)

#-------------------------------------------------------------------------------
#Viewing and removing completely empty rows

empty_rows <- HPV_genotyping_data %>%
  filter(if_all(everything(), is.na))

# view count of empty_rows then
# remove from main data frame #for viewing purposes (optional)

#cat("Dataset consist of", nrow(empty_rows), " empty rows")

HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(!if_all(everything(), is.na)) 
#-------------------------------------------------------------------------------
# Viewing and removing duplicated rows

duplicated_rows <- HPV_genotyping_data %>%
  filter(duplicated(.))

# view count of duplicated_rows then
# remove from main data frame #for viewing purposes (optional)

#cat("Dataset consist of", nrow(duplicated_rows), " duplicated rows")

HPV_genotyping_data <- HPV_genotyping_data %>%
  distinct()

# Check if duplicates remain

if (any(duplicated(HPV_genotyping_data))) {
  cat("There are still duplicated rows in the dataframe.\n")
} else {
  cat("All duplicated rows have been removed.\n")
}

#-------------------------------------------------------------------------------
# Keep only records with dates up to and including 2022

HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(large_blood_date <= as.Date("2022-12-31"))

#-------------------------------------------------------------------------------
# Clean column names first to ensure consistency using clean_names()

HPV_genotyping_data <-HPV_genotyping_data %>% 
  clean_names()

#-------------------------------------------------------------------------------
# Renaming HPV type columns originally labeled as x_(number)
# to a clearer format like HPV16, HPV18, etc., for readability.
# and add zero to single-digit HPV types

names(HPV_genotyping_data) <- names(HPV_genotyping_data) %>%
  gsub("^x(\\d+)$", "HPV\\1", .) %>%   
  gsub("^HPV(\\d)$", "HPV0\\1", .)      



#-------------------------------------------------------------------------------
# Sort dataset first by study_id, then by visit date
# This helps us see individuals with multiple visits and the dates they came in

HPV_genotyping_data <- HPV_genotyping_data %>%
  arrange(study_id, large_blood_date)

#-------------------------------------------------------------------------------
# identifying individuals with >3 visits and removing them

more_than_3_visits <- HPV_genotyping_data %>%
  distinct(study_id, large_blood_date) %>%
  count(study_id, name = "num_visits") %>%
  filter(num_visits > 3)

HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(!study_id %in% more_than_3_visits$study_id)

#-------------------------------------------------------------------------------

# Check for unique values in each column
# except study_id, visit_code, and large_blood_date

cols_to_check <- setdiff(names(HPV_genotyping_data),
                         c("study_id", "large_blood_date", "visit_code"))

for (col in cols_to_check) {
  cat("\n------", col, "------\n")
  print(unique(HPV_genotyping_data[[col]]))
}


#-------------------------------------------------------------------------------
# Following are changes based on unique (from above)
# In comments col :Change "Blood contaminatated" to "Blood contaminated"

HPV_genotyping_data$comment <- gsub("Blood contaminatated",
                                    "Blood contaminated",
                                    HPV_genotyping_data$comment)

# Check if the incorrect spelling still exists
incorrect_rows <- HPV_genotyping_data %>%
  filter(grepl("Blood contaminatated", comment)) %>%
  select(comment)

if (nrow(incorrect_rows) == 0) {
  cat("No rows with 'Blood contaminatated' found — replacement successful!")
} else {
  cat("Rows with 'Blood contaminatated' still present:")
  print(incorrect_rows)
}

#-------------------------------------------------------------------------------
# Replace string "NA" with actual NA values across character and 
# factor columns to ensure proper missing data handling 
# thus Standardize missing values

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
  mutate(across(where(is.factor), ~ na_if(., "NA")))

#-------------------------------------------------------------------------------
# view data type of each column 

str(HPV_genotyping_data)

#-------------------------------------------------------------------------------
# Data type conversion

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(
    large_blood_date = as.Date(large_blood_date),
    across(c(nugent_gardnerella,nugent_mobiluncus ,nugent_lactobacillus ,
             nugent_total_score), as.numeric))
    
#-------------------------------------------------------------------------------
# Column 1 check - inspecting pattern of study_id
# Flag invalid study_id values (wrong format or missing)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(
    study_id_valid = if_else(
      !is.na(study_id) & grepl("^\\d{3}-\\d{2}-\\d{4}-\\d{3,4}$", study_id),
      TRUE, FALSE))

# Extract invalid records for inspection
invalid_study_ids <- HPV_genotyping_data %>%
  filter(study_id_valid == FALSE) %>%
  mutate(flag_reason = "Invalid or missing study_id") %>%
  select(study_id,flag_reason,study_id_valid)

# Keep only valid IDs for main analyses and remove temporary flag column
HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(study_id_valid == TRUE)%>%
  select(-study_id_valid)

#-------------------------------------------------------------------------------
#column 2 check - inspecting pattern
# Flag invalid study_id values (wrong format) 

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(
    invalid_code = if_else(grepl("^\\d{1}-\\d{3}-\\d{1}$", visit_code),
      TRUE, FALSE))

# Extract invalid records for inspection
invalid_visit_codes <- HPV_genotyping_data %>%
  filter(invalid_code == FALSE) %>%
  mutate(flag_reason = "Invalid visit code") %>%
  select(study_id,flag_reason,invalid_code)

# Keep only valid visit codes for main analyses and remove temporary flag column
HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(invalid_code == TRUE)%>%
  select(-invalid_code)


#-------------------------------------------------------------------------------
#column 7 check -
#check how many rows have incorrect Nugent total scores
#a mismatch is where the nugent_total_score doesn't equal the sum of its components

mismatch_count <- sum(!is.na(HPV_genotyping_data$nugent_total_score) &
    HPV_genotyping_data$nugent_total_score != (HPV_genotyping_data$nugent_gardnerella +
                                                 HPV_genotyping_data$nugent_mobiluncus +
                                                 HPV_genotyping_data$nugent_lactobacillus))

#corrects miscalculated rows
#only updates rows where the total is incorrect

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(nugent_total_score = if_else(
    !is.na(nugent_total_score) &
      nugent_total_score != (nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus),
    nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus,
    nugent_total_score))

cat("Nugent total scores corrected for", mismatch_count, "rows.")

#-------------------------------------------------------------------------------
#column 8 check -
# Count mismatched interpretations

interpret_mismatch_count <- HPV_genotyping_data %>%
  filter(!is.na(nugent_total_score) & !is.na(nugent_interpretation) &
           (
             (nugent_total_score < 3 & nugent_interpretation != "NO BV") |
               (nugent_total_score >= 4 & nugent_total_score <= 6 & nugent_interpretation != "INTERMEDIATE") |
               (nugent_total_score >= 7 & nugent_interpretation != "BV")
           )
  ) %>% nrow()

# corrects mismatched interpretations

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(nugent_interpretation = case_when(
    nugent_total_score < 3 ~ "NO BV",
    nugent_total_score >= 4 & nugent_total_score <= 6 ~ "INTERMEDIATE",
    nugent_total_score >= 7 ~ "BV",
    TRUE ~ nugent_interpretation  
  ))

cat("Nugent interpretations corrected for", interpret_mismatch_count, "rows.")

#-------------------------------------------------------------------------------
# column check 9-12 - converting  "DETECTED" and "NOT DETECTED" text to binary
# to keep consistent reporting format for both STI and HPV results
# where "DETECTED" ~ 1 (yes there is an sti)
# where "NOT DETECTED" ~ 0 (no there is no sti)

sti_cols <- grep("^sti_", names(HPV_genotyping_data), value= TRUE)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(across(all_of(sti_cols),
                ~ case_when(. == "DETECTED" ~ 1,
                            . == "NOT DETECTED" ~ 0,
                            TRUE ~ NA_real_ )))
#-------------------------------------------------------------------------------
# Replicates - concordant and discordant individuals -
# Creating a replicate data

replicate_data <- HPV_genotyping_data %>%
  group_by(study_id, large_blood_date) %>%
  filter(n() >1) %>%    #returns individuals with more than 1 sample
  ungroup()

# Classify each replicate visit as concordant or discordant

hpv_type_col <- grep("HPV\\d+", names(HPV_genotyping_data), value = TRUE)

replicate_classification <- replicate_data %>%
  group_by(study_id, large_blood_date) %>%
  summarise(is_concordant = n_distinct(across(all_of(hpv_type_col))) == 1 , .groups = "drop")

# data frames to work on separately - concordant_visits and discordant visits 

concordant_visits <- replicate_data %>%
  semi_join(replicate_classification %>% 
              filter(is_concordant),
            by = c("study_id", "large_blood_date"))

discordant_visits <- replicate_data %>%
  semi_join(replicate_classification %>% 
              filter(!is_concordant),
            by = c("study_id", "large_blood_date"))


# vectors to store cols from replicate data 

nugent_col <- grep("nugent_\\D", names(replicate_data), value = TRUE)

# Combine all relevant variable groups (Nugent, STI, and HPV) into one vector.
# This allows simultaneous NA-counting across all key variables instead of checking each group separately.

all_info_cols <- c(nugent_col, sti_cols, hpv_type_col)

# for concordant individuals remove the other "same" case :differences are with
# BV and STI status keep rows with fewest NA values 

cleaned_concordant_replicates <- concordant_visits %>%
  group_by(study_id, large_blood_date) %>%
  # Count NA per row across relevant columns
  mutate(na_count = rowSums(is.na(across(all_of(all_info_cols))))) %>%
  # Keep the row with the fewest NAs per replicate group
  slice_min(order_by = na_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-na_count)


# For discordant individuals (same study_id and large_blood_date), 
# retain any positive HPV results (1) across all rows, 
# and preserve the first available non-NA result for Nugent and STI row 

collapsed_discordant_data <- discordant_visits %>%
  group_by(study_id, visit_code, large_blood_date) %>%
  summarise(
    
    # Collapse Nugent and STI columns (first non-NA)
    across(all_of(c(nugent_col, sti_cols)), ~ {
      valid <- .[!is.na(.)]
      if (length(valid) > 0) valid[1] else NA
    }, .names = "{.col}"),
    
    # Collapse HPV columns (max, treating NA safely)
    across(all_of(hpv_type_col), ~ {
             if (all(is.na(.))) NA_integer_ else max(., na.rm = TRUE)}, .names = "{.col}"),
    
    # report: retain "result" if present, else NA
    report = {valid <- report[report == "RESULT"]
      if (length(valid) > 0) valid[1] else NA_character_},
    
    # comment: ONLY keep NA (drop all non-NA)
    comment = {if (all(is.na(comment))) NA_character_ else NA_character_ },
    .groups = "drop"
    )


# anti_join replicates then bind_rows the fixed non-replicates 

HPV_genotyping_data <- anti_join(HPV_genotyping_data, replicate_data)

HPV_genotyping_data <- bind_rows(HPV_genotyping_data,cleaned_concordant_replicates)

HPV_genotyping_data <- bind_rows(HPV_genotyping_data,collapsed_discordant_data)

# arrange again since bind_rows added rows from the last observation 

HPV_genotyping_data <- HPV_genotyping_data %>%
  arrange(study_id, large_blood_date)
  

#-------------------------------------------------------------------------------
#Column 8 is categorical where we encode set numbers 0,1,2
# 0 = "NO BV"
# 1 = "BV"
# 2 = "INTERMEDIATE"

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(nugent_interpretation = case_when
         (is.na(nugent_interpretation) ~ NA_real_,
           nugent_interpretation == "NO BV" ~ 0,
           nugent_interpretation == "INTERMEDIATE" ~ 2,
           nugent_interpretation == "BV" ~ 1))
#-------------------------------------------------------------------------------
# new added columns
# 1. HPV_sum (Total no of HPV an individual is infected with)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(HPV_sum = rowSums(select(., all_of(hpv_type_col)), na.rm = FALSE))

#------------------------------------------------------------------------------
# 2. HR-HPV_sum (Total no of High-risk HPV an individual is infected with )

hr_hpv_types <- c("HPV16", "HPV18", "HPV26", "HPV31", "HPV33", "HPV35",
                  "HPV39", "HPV45", "HPV51", "HPV52", "HPV53", "HPV56",
                  "HPV58", "HPV59", "HPV66", "HPV68", "HPV73", "HPV82")

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(HR_HPV_sum = rowSums(select(., all_of(hr_hpv_types)), na.rm = TRUE))

#-------------------------------------------------------------------------------
# 3.LR-HPV_sum (Total no of Low-risk HPV an individual is infected with)

lr_hpv_types <- c("HPV06", "HPV11", "HPV40", "HPV42", "HPV43",
                  "HPV54", "HPV61", "HPV67", "HPV69",
                  "HPV70", "HPV71", "HPV72", "HPV84")

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(LR_HPV_sum = rowSums(select(., all_of(lr_hpv_types)), na.rm = TRUE))

#-------------------------------------------------------------------------------
# 4.STI_pcr_count(Total no of number of infected sti to see relation between the
# number and status of HPV type(what does it mean?))

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(sti_pcr_count = rowSums(select(., all_of(sti_cols)), na.rm = TRUE))

#-------------------------------------------------------------------------------
# 5.coinfection (with other STI)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(coinfection = case_when(
    !is.na(HPV_sum) & !is.na(sti_pcr_count) & HPV_sum > 0 & sti_pcr_count > 0 ~ 1,
    !is.na(HPV_sum) & !is.na(sti_pcr_count) & (HPV_sum == 0 | sti_pcr_count == 0) ~ 0,
    TRUE ~ NA_real_ ))

#-------------------------------------------------------------------------------
# 6.HPV group (classified to understand multiple infection in individuals)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(HPV_group = case_when(
    is.na(HPV_sum) ~ NA_character_,
    HPV_sum == 0 ~ "No infection",
    HPV_sum == 1 ~ "Single infection",
    HPV_sum >= 2 ~ "Multiple infection"))

#-------------------------------------------------------------------------------
# 7.HPV
# (if an individual has HPV =1
# if an individual does not have HPV)

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(HPV_positive = case_when(
      rowSums(select(., all_of(hpv_type_col)) == 1, na.rm = TRUE) > 0 ~ 1,
      rowSums(!is.na(select(., all_of(hpv_type_col)))) == 0 ~ NA_real_,
      TRUE ~ 0))
#-------------------------------------------------------------------------------
# 8. HR_HPV_positive if individual has any high-risk HPV type = 1
#    LR_HPV_positive if individual has any low-risk HPV type = 1
#    otherwise : no = 0
#    NA remains NA

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(
    HR_HPV_positive = case_when(
      rowSums(!is.na(select(., all_of(hr_hpv_types)))) == 0 ~ NA_real_,
      rowSums(select(., all_of(hr_hpv_types)) == 1, na.rm = TRUE) > 0 ~ 1,
      TRUE ~ 0
    ),
    
    LR_HPV_positive = case_when(
      rowSums(!is.na(select(., all_of(lr_hpv_types)))) == 0 ~ NA_real_,
      rowSums(select(., all_of(lr_hpv_types)) == 1, na.rm = TRUE) > 0 ~ 1,
      TRUE ~ 0
    ))

#-------------------------------------------------------------------------------
# 9.mixed_infection
# indicates whether an individual is infected with:
# both HR-HPV and LR-HPV = 1
# if only infected with the one = 0 (no mixed infection)
# NA remains NA 

HPV_genotyping_data <- HPV_genotyping_data %>%
  mutate(mixed_risk_infection = case_when(
    is.na(HR_HPV_sum) | is.na(LR_HPV_sum) ~ NA_real_,
    HR_HPV_sum > 0 & LR_HPV_sum > 0 ~ 1,
    TRUE ~ 0))


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# clean data set!
 
HPV_genotyping_data <- HPV_genotyping_data %>%
  select(
    study_id,
    visit_code,
    large_blood_date,
    nugent_gardnerella,
    nugent_mobiluncus,
    nugent_lactobacillus,
    nugent_total_score,
    nugent_interpretation,
    sti_mycoplasma_genitalium_pcr,
    sti_trichomonas_vaginalis_pcr,
    sti_neisseria_gonorrhoeae_pcr,
    sti_chlamydia_trachomatis_pcr,
    sti_pcr_count,
    HPV_positive,
    HPV16, HPV18, HPV26, HPV31, HPV33, HPV35, HPV39, HPV45,
    HPV51, HPV52, HPV53, HPV56, HPV58, HPV59, HPV66, HPV68,
    HPV73, HPV82,
    HR_HPV_positive,
    HR_HPV_sum,
    HPV06, HPV11, HPV40, HPV42, HPV43, HPV54,
    HPV61, HPV67, HPV69, HPV70, HPV71, HPV72, HPV84,
    LR_HPV_positive,
    LR_HPV_sum,
    HPV_sum,
    mixed_risk_infection,
    HPV_group,
    coinfection,
    report)

# Table with no NA values across HPV columns 
# Final data set will have individuals with complete genotyping data
# if_all keep only rows where all HPV columns are filled
# mutate NEW column called visit_num 

HPV_genotyping_complete <- HPV_genotyping_data %>%
  filter(if_all(all_of(hpv_type_col), ~ !is.na(.)))%>%
  group_by(study_id) %>%
  mutate(visit_num = row_number()) %>%  
  ungroup()

# Final data set here includes all individuals even those with NA for genotying data 
# NOT TO BE USED FOR ANALYSIS but for reporting missing data
HPV_genotyping_data_final_clean <- HPV_genotyping_data %>%
  group_by(study_id) %>%
  mutate(visit_num = row_number()) %>%  
  ungroup()

# Save cleaned data set as csv
# this dataset contains all individuals NA hpv tests and those with results  
write_csv(HPV_genotyping_data, "HPV_genotyping_data_final_clean.csv")

# this dataset is what we use for all analysis -- does not contain NA VALUES 
# for HPV represents true tested population 
write_csv(HPV_genotyping_complete, "HPV_genotyping_complete.csv")
#-------------------------------------------------------------------------------