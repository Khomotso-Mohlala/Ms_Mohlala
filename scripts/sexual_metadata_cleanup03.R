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
library(writexl)     #Write Excel file
library(lubridate)   #Date formatting and manipulation

#-------------------------------------------------------------------------------
#set working directory to folder that contains excel files

WKDIR <- "C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts"
setwd(WKDIR)

#-------------------------------------------------------------------------------
#load data into R using excel file name and view table  

filename <- "03. FRESH_HPV_sexual_behavior_metadata_04.06.2025.xlsx"

# Check if file exists before attempting to read it

if (file.exists(filename)) {
  data_03  <- read_excel(filename, sheet = 'Data')
  View(data_03)
} else {
  stop("File not found: ", filename)
}

#-------------------------------------------------------------------------------  
# Keep only key variables for cleaning 
# Rename column names

Sexual_behaviour_metadata <- data_03 %>%
  select( 
    study_id = subjid,
    visit_id = visitid,
    visit_date = vstdt,
    visit_code,
    hiv_risk_date_3month = visit_date_hiv_risk_3month,
    sex_episodes_7d = numsexepisodes7day,
    sex_partners_7d = numsexpartners7days,
    sex_episodes_30d = numsexepisodes30day,
    condoms_used_30d = condomsused30days,
    sex_partners_30d = numsexpartners30day,
    oral_sex_7d = oral_sex_last_7_days_yn,
    penile_vaginal_sex_7d = penile_vaginal_last_7_days_yn,
    anal_sex_7d = anal_last_7_days_yn,
    family_planning_type = familyplanningtype,
    sti_symptoms_now = stisymptomsnow,
    date_lower,
    date_upper,
    age_at_baseline = age_baseline,
    age_first_sex_baseline = age_first_sex_baseline_demo,
    number_sexual_partners,
    live_with_partner_baseline)


#-------------------------------------------------------------------------------
#Viewing and removing completely empty rows 

empty_rows <- Sexual_behaviour_metadata %>% 
  filter(if_all(everything(), is.na))

# view count of empty_rows then
# Drop empty_rows 

# cat("Dataset consist of", nrow(empty_rows), " empty rows")

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  filter(!if_all(everything(), is.na)) 

#-------------------------------------------------------------------------------
# Viewing and removing duplicated rows 

duplicated_rows <- Sexual_behaviour_metadata %>%
  filter(duplicated(.))

# view count of duplicated_rows then
# Drop duplicated rows 

#cat("Dataset consist of", nrow(duplicated_rows), " duplicated rows")

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  distinct()

# Check if duplicates remain

if (any(duplicated(Sexual_behaviour_metadata))) {
  cat("There are still duplicated rows in the dataframe.\n")
} else {
  cat("All duplicated rows have been removed.\n")
}

#-------------------------------------------------------------------------------
# dataframe with individuals from before and including   2022

Sexual_behaviour_metadata <- Sexual_behaviour_metadata  %>%
  filter(visit_date <= as.Date("2022-12-31"))

#-------------------------------------------------------------------------------
# Clean column names first to ensure consistency using clean_names() 

Sexual_behaviour_metadata <-Sexual_behaviour_metadata %>% 
  clean_names()
  
#-------------------------------------------------------------------------------
# Sort dataset first by study_id, then by visit date
# This helps us see individuals with multiple visits and the dates they came in 

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  arrange(study_id, visit_date)  

#-------------------------------------------------------------------------------
# Replace string "NA" with actual NA values across character and 
# factor columns to ensure proper missing data handling 
# thus Standardize missing values

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
  mutate(across(where(is.factor), ~ na_if(., "NA")))


#-------------------------------------------------------------------------------
# view data type of each column 

str(Sexual_behaviour_metadata)

#-------------------------------------------------------------------------------
# data type conversion 

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(across(c(visit_date, hiv_risk_date_3month, date_upper, date_lower), as.Date),
         across(c(sex_episodes_7d, sex_partners_7d, sex_episodes_30d, sex_partners_30d,
                   condoms_used_30d, number_sexual_partners,live_with_partner_baseline,
                  oral_sex_7d,penile_vaginal_sex_7d,anal_sex_7d,sti_symptoms_now,
                  live_with_partner_baseline,family_planning_type),
                 as.numeric))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# column 1 check -

# Flag invalid study_id values (wrong format or missing)
Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(
    study_id_valid = if_else(
      !is.na(study_id) & grepl("^\\d{3}-\\d{2}-\\d{4}-\\d{3,4}$", study_id),
      TRUE, FALSE))

# Extract invalid records for inspection
invalid_study_ids <- Sexual_behaviour_metadata %>%
  filter(study_id_valid == FALSE) %>%
  mutate(flag_reason = "Invalid or missing study_id") %>%
  select(study_id,flag_reason,study_id_valid)


# Keep only valid IDs for main analyses and remove temporary flag column
Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  filter(study_id_valid == TRUE)%>%
select(-study_id_valid)


#-------------------------------------------------------------------------------
# column 2 and 4 check -
# inspecting similarity in visit_id and visit_code 

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(
    same_codes = case_when(
      is.na(visit_id) | is.na(visit_code) ~ 2,       
      visit_id == visit_code ~ 1,                   
      TRUE ~ 0  ))                                     

# mismatch
mismatch <- Sexual_behaviour_metadata %>%
  filter(same_codes == 0)

# at least one NA 
one_na_value <- Sexual_behaviour_metadata %>%
  filter(same_codes == 2)

# match
match <- Sexual_behaviour_metadata %>% 
  filter(same_codes == 1)

# mismatch accounts for 45/902 obs. accounting for 4.99% therefore keep visit_id 
# and since there is incomplete observations in visit-code exclude from final dataset

#-------------------------------------------------------------------------------
# column 21 check - 

unique(Sexual_behaviour_metadata$live_with_partner_baseline)

# Encode live_with_partner_baseline: 1 = yes, 2 = no make 1/0

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(live_with_partner_baseline = case_when(
    live_with_partner_baseline == 1 ~ 1,
    live_with_partner_baseline == 2 ~ 0,
    TRUE ~ NA_real_))
  
#-------------------------------------------------------------------------------
# column age check -
# change age data type from chr to num

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(across(c(age_at_baseline, age_first_sex_baseline ),
                as.numeric))

# To inspect age ranges thus -- accepted range 18-25 years(age_at_baseline)

summary(Sexual_behaviour_metadata$age_at_baseline)
unique(Sexual_behaviour_metadata$age_at_baseline)

summary(Sexual_behaviour_metadata$age_first_sex_baseline)
unique(Sexual_behaviour_metadata$age_first_sex_baseline)


#-------------------------------------------------------------------------------
# column dates check-
# To observe the relationships between dates 

# showing time gaps in days
Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(range_lower_to_upper = as.numeric(date_upper - date_lower),)
  
# Extract records where range > 28 days for inspection
# As expected, the length between first visit and follow-up is approximately one month (≈28 days),
# so any value exceeding 28 days is flagged for further review.

invalid_date_ranges <- Sexual_behaviour_metadata %>%
  filter(!is.na(range_lower_to_upper) & range_lower_to_upper > 28) %>%
  mutate(flag_reason = paste0("Date range exceeds 28 days (", range_lower_to_upper, " days)")) %>%
  select(study_id, flag_reason, range_lower_to_upper)


# Keep only valid date records in main dataset
# Remove flagged entries and drop temporary column.
Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  filter(is.na(range_lower_to_upper) | range_lower_to_upper <= 28) %>%
  select(-range_lower_to_upper)
#-------------------------------------------------------------------------------
# new added columns
# 1.years_since_first_sex (an estimate of how long (in years) each individual has been sexually active)

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(
    years_since_first_sex = as.numeric(age_at_baseline) - as.numeric(age_first_sex_baseline))

#-------------------------------------------------------------------------------
# creating  a variable that reports if participant had 
# any type of sex in the last 7 days.
# 1= yes 0= no

Sexual_behaviour_metadata <- Sexual_behaviour_metadata %>%
  mutate(
    had_sex_7d = case_when(
      oral_sex_7d == 1 | penile_vaginal_sex_7d == 1 | anal_sex_7d == 1 ~ 1,
      oral_sex_7d == 0 & penile_vaginal_sex_7d == 0 & anal_sex_7d == 0 ~ 0,
      TRUE ~ NA_real_  )) 
#-------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------  
#clean data set!

Sexual_behaviour_metadata <- Sexual_behaviour_metadata%>%
  select(
    study_id,
    visit_id,
    visit_date,
    sex_episodes_7d,
    sex_partners_7d,
    sex_episodes_30d,
    sex_partners_30d,
    condoms_used_30d,
    number_sexual_partners,
    oral_sex_7d,
    penile_vaginal_sex_7d,
    anal_sex_7d,
    had_sex_7d,
    family_planning_type,
    sti_symptoms_now,
    age_at_baseline,
    age_first_sex_baseline,
    years_since_first_sex,
    live_with_partner_baseline)

# Save cleaned data set as CSV

write_csv(Sexual_behaviour_metadata, "Sexual_behaviour_metadata_final_clean.csv")

# Export inspection log for later review

Sexual_behaviour_metadata_inspection <- bind_rows(invalid_date_ranges,invalid_study_ids)

#-------------------------------------------------------------------------------







 













