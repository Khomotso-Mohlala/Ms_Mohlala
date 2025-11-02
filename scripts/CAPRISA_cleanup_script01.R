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

library(readxl)      # Read Excel files
library(tidyverse)   # Tidy data tools 
library(janitor)     # Clean data frames and column names
library(writexl)     # Write Excel files
library(lubridate)   # Date formatting and manipulation

# -------------------------------------------------------------------------------
#set working directory to folder that contains excel files

WKDIR <- "C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts"
setwd(WKDIR)

#-------------------------------------------------------------------------------
#load data into R from Excel file  

filename <- "01_First_delivery_at_caprisa_all_samples.xlsx"

# Check if file exists before attempting to read it

if (file.exists(filename)) {
  data_01 <- read_excel(filename, sheet = "Sheet1")
  View(data_01)
} else {
  stop("File not found: ", filename)
}

#-------------------------------------------------------------------------------
# Keep only key variables for cleaning

CAPRISA_data_01 <- data_01 %>%
  select(
    study_id,
    visit_code,
    large_blood_date,
    study,
    suuba_date,
    sample,
    location,
    num,
    nugent_gardnerella,
    nugent_mobiluncus,
    nugent_lactobacillus,
    nugent_total_score,
    nugent_interpretation,
    sti_mycoplasma_genitalium_pcr,
    sti_trichomonas_vaginalis_pcr,
    sti_neisseria_gonorrhoeae_pcr,
    sti_chlamydia_trachomatis_pcr)

View(CAPRISA_data_01)

#-------------------------------------------------------------------------------
#Viewing and removing completely empty rows 

empty_rows <- CAPRISA_data_01 %>% 
  filter(if_all(everything(), is.na))

#View(empty_rows)  #for viewing purposes (optional) 

cat("Empty rows removed:", nrow(empty_rows))

CAPRISA_data_01 <- anti_join(CAPRISA_data_01, empty_rows)

#-------------------------------------------------------------------------------
# Viewing and removing duplicated rows 

duplicated_rows <- CAPRISA_data_01 %>%
  filter(duplicated(.))

# View(duplicated_rows) #for viewing purposes (optional)

cat("Duplicates removed:", nrow(duplicated_rows))

CAPRISA_data_01 <- anti_join(CAPRISA_data_01, duplicated_rows)

#-------------------------------------------------------------------------------
# Clean column names first to ensure consistency using clean_names() 

CAPRISA_data_01 <- CAPRISA_data_01 %>% 
  clean_names()

#-------------------------------------------------------------------------------
# Convert string "NA" to proper NA to avoid hidden missing data

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(across(where(is.character), ~ na_if(., "NA"))) %>%
  mutate(across(where(is.factor), ~ na_if(., "NA")))

#-------------------------------------------------------------------------------
# Check data type of all columns

str(CAPRISA_data_01)

#-------------------------------------------------------------------------------
# Data type conversion

CAPRISA_data_01 <- CAPRISA_data_01 %>% 
  mutate(
    nugent_gardnerella = as.numeric(nugent_gardnerella),
    nugent_mobiluncus = as.numeric(nugent_mobiluncus),
    nugent_lactobacillus = as.numeric(nugent_lactobacillus),
    nugent_total_score = as.numeric(nugent_total_score))
#-------------------------------------------------------------------------------
# Sort data set first by study_id, then by visit date
# This helps us see individuals with multiple visits and the dates they came in 

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  arrange(study_id, large_blood_date)

#-------------------------------------------------------------------------------
# identifying individuals with >3 visits and removing them

more_than_3_visits <- CAPRISA_data_01 %>%
  distinct(study_id, large_blood_date) %>%
  count(study_id, name = "num_visits") %>%
  filter(num_visits > 3)

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  filter(!study_id %in% more_than_3_visits$study_id)

#-------------------------------------------------------------------------------
# Column 1 check - inspecting pattern of study_id

invalid_ids <- CAPRISA_data_01 %>%
  filter(!grepl("^\\d{3}-\\d{2}-\\d{4}-\\d{3,4}$", study_id))

print(invalid_ids$study_id)

#-------------------------------------------------------------------------------
#column 2 check - inspecting pattern 

invalid_code <- CAPRISA_data_01 %>% 
  filter(!grepl("^\\d{1}-\\d{3}-\\d{1}$", visit_code))

print(invalid_code$visit_code)

#-------------------------------------------------------------------------------
#column 3 check - standardizing date to YY-MM-DATE

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(large_blood_date = mdy(large_blood_date))

#-------------------------------------------------------------------------------
# Keep only records with dates up to and including 2022

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  filter(large_blood_date <= as.Date("2022-12-31"))
# column 4 check -

unique(CAPRISA_data_01$study)
#-------------------------------------------------------------------------------
# column 5 check -
# suuba_date is indicative of a date format stored in excel where it calculates 
# the number of days since January 1, 1900 
# convert this to a readable date format

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(suuba_date = as.Date(as.numeric(suuba_date), origin = "1899-12-30"))

# checking for similarity match between suuba_date and large_blood_date 

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(same_dates = large_blood_date == suuba_date) 

#Print where case if false where large_blood_date != suuba_date

CAPRISA_data_01 %>%
  filter(same_dates == FALSE)

View(CAPRISA_data_01)
#-------------------------------------------------------------------------------
# column 6 check -

unique(CAPRISA_data_01$sample)
 
CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(sample = gsub("CVL_PELLET", "CVL_Pellet", sample))

#-------------------------------------------------------------------------------
# column 8 check -

unique(CAPRISA_data_01$location)
#-------------------------------------------------------------------------------
# column 11 check -

unique(CAPRISA_data_01$num)
#-------------------------------------------------------------------------------
# column nugent check -
# check how many rows have incorrect Nugent total scores
# a mismatch is where the nugent_total_score doesn't equal the sum of its 
# components

mismatch_count <- sum(
  !is.na(CAPRISA_data_01$nugent_total_score) &
    CAPRISA_data_01$nugent_total_score != (CAPRISA_data_01$nugent_gardnerella +
                                                 CAPRISA_data_01$nugent_mobiluncus +
                                                 CAPRISA_data_01$nugent_lactobacillus))

#corrects miscalculated rows 
#only updates rows where the total is incorrect
CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(nugent_total_score = if_else(
    !is.na(nugent_total_score) &
      nugent_total_score != (nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus),
    nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus,
    nugent_total_score
  ))

cat("Nugent total scores corrected for", mismatch_count, "rows.")

#-------------------------------------------------------------------------------
# column nugent_interpretation check -
# Count mismatched interpretations 

interpret_mismatch_count <- CAPRISA_data_01 %>%
  filter(!is.na(nugent_total_score) & !is.na(nugent_interpretation) &
           (
             (nugent_total_score < 3 & nugent_interpretation != "NO BV") |
               (nugent_total_score >= 4 & nugent_total_score <= 6 & nugent_interpretation != "INTERMEDIATE") |
               (nugent_total_score >= 7 & nugent_interpretation != "BV")
           )
  ) %>% nrow()

# corrects mismatched interpretations

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(nugent_interpretation = case_when(
    nugent_total_score < 3 ~ "NO BV",
    nugent_total_score >= 4 & nugent_total_score <= 6 ~ "INTERMEDIATE",
    nugent_total_score >= 7 ~ "BV",
    TRUE ~ nugent_interpretation  
  ))

cat("Nugent interpretations corrected for", interpret_mismatch_count, "rows.")
#-------------------------------------------------------------------------------
# Column nugent_interpretation is categorical where we encode set numbers 0,1,2,3,4
# 0 = "NO BV" 
# 1 = "BV"
# 2 = "INTERMEDIATE"
# 4 = NA

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(nugent_interpretation = case_when 
         (is.na(nugent_interpretation) ~ 4,
           nugent_interpretation == "NO BV" ~ 0,
           nugent_interpretation == "INTERMEDIATE" ~ 2,
           nugent_interpretation == "BV" ~ 1))

#-------------------------------------------------------------------------------
# column check 9-12 - converting  "DETECTED" and "NOT DETECTED" text to binary
# to keep consistent reporting format for both STI and HPV results
# where "DETECTED" ~ 1 (yes there is an sti)
# where "NOT DETECTED" ~ 0 (no there is no sti)

sti_cols <- c(
  "sti_mycoplasma_genitalium_pcr", 
  "sti_trichomonas_vaginalis_pcr",
  "sti_neisseria_gonorrhoeae_pcr", 
  "sti_chlamydia_trachomatis_pcr")

CAPRISA_data_01 <- CAPRISA_data_01 %>%
  mutate(across(all_of(sti_cols),
                ~ case_when(. == "DETECTED" ~ 1,
                            . == "NOT DETECTED" ~ 0,
                            TRUE ~ NA_real_)))
                
View(CAPRISA_data_01)
#-------------------------------------------------------------------------------
# clean data set!

CAPRISA_data_01 <- CAPRISA_data_01 %>% 
  select(
    study_id,
    visit_code,
    large_blood_date,
    study,
    suuba_date,
    sample,
    location,
    nugent_gardnerella,
    nugent_mobiluncus,
    nugent_lactobacillus,
    nugent_total_score,
    nugent_interpretation,
    sti_mycoplasma_genitalium_pcr,
    sti_trichomonas_vaginalis_pcr,
    sti_neisseria_gonorrhoeae_pcr,
    sti_chlamydia_trachomatis_pcr)


# Save cleaned data set as CSV

write_csv(CAPRISA_data_01, "CAPRISA_data_01.csv")

#-------------------------------------------------------------------------------

























