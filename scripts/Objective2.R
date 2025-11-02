#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Objective 2 :To assess the prevalence and patterns of co-persistent HPV 
#              infections.

# Definition:
#                           
# 1. Co-persistent HPV infection refers to the persistence of two or more of the 
# same HPV genotypes across three consecutive visits. 

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Libraries to load for OBJECTIVE 2
#-------------------------------------------------------------------------------
# Set up of libraries needed
# To download type install.packages("") name of library as below in "quotes"

library(tidyverse)   # collection of tools for data manipulation and visualization
library(pheatmap)    # for heatmaps
library(knitr)       # for knitting tables/reports
library(kableExtra)  # for publication-ready tables 
library(writexl)     # for exporting to Excel
#-------------------------------------------------------------------------------
#set path to working directory

WKDIR <- "C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts"
setwd(WKDIR)

# Read in csv
# dataset excludes individuals with missing HPV results (including full genotyping info)
df_04 <- read_csv("HPV_genotyping_complete.csv")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# create vectors to store HPV genotype columns 
# first line is all columns (all HPV types in data frame) 
# second line is high risk types only 
# third line is low risk types only

HPV_cols <- grep("HPV\\d+", names(df_04), value= TRUE)

hr_hpv_types <- c("HPV16", "HPV18", "HPV26", "HPV31", "HPV33", "HPV35",
                  "HPV39", "HPV45", "HPV51", "HPV52", "HPV53", "HPV56",
                  "HPV58", "HPV59", "HPV66", "HPV68", "HPV73", "HPV82")

lr_hpv_types <- c("HPV06", "HPV11", "HPV40", "HPV42", "HPV43",
                  "HPV54", "HPV61", "HPV67", "HPV69",
                  "HPV70", "HPV71", "HPV72", "HPV84")

#------------------------------------------------------------------------------#
# Prepare data frame to store individuals with 3 visits 

complete_visits <- df_04 %>%
  group_by(study_id) %>%
  filter(n_distinct(visit_num) == 3, all(c(1, 2, 3) %in% visit_num)) %>%
  ungroup() %>%
  select(study_id, all_of(HPV_cols), HPV_sum, visit_num)

# save csv for future use
write_csv(complete_visits, file = "complete_consecutive_visits.csv")

# Use collapse to determine which genotypes persist across all 3 visits
# thus creating one entry per individual with persisting HPV genotypes 

persistent_hpv <- complete_visits %>%
  group_by(study_id) %>%
  summarise(
    across(all_of(HPV_cols), ~ as.integer(all(.x == 1))),
    .groups = "drop")

# save csv for future use 
write_csv(persistent_hpv, file = "persistening_hpv_patterns.csv")

# Determine each individual's multiple persistent HPV genotype infection profile 
# Summarize the frequency, percentage, and risk category (HR, LR, or Mixed) of these patterns

patterns <- persistent_hpv %>%
  # keep only rows with 2+ persistent types
  mutate(n_pos = rowSums(select(., all_of(HPV_cols)), na.rm = TRUE)) %>%
  filter(n_pos >= 2) %>%
  
  # build the pattern string from the kept HPV columns
  mutate(
    persistent_pattern = apply(select(., all_of(HPV_cols)), 1, function(x) {
      pos <- names(x)[x == 1]
      paste(pos, collapse = " + ")})) %>%
  select(study_id, persistent_pattern) %>%
  count(persistent_pattern, name = "Pattern_frequency") %>%
  mutate(percent_persistent = round(100 * Pattern_frequency / sum(Pattern_frequency), 2)) %>%
  
  # risk classification
  mutate(
    types_list = str_split(persistent_pattern, " \\+ ") %>% map(~ str_trim(.x)),
    risk_category_persistent = map_chr(types_list, function(types) {
      if (length(types) == 0) "None"
      else if (all(types %in% hr_hpv_types)) "HR only"
      else if (all(types %in% lr_hpv_types)) "LR only"
      else "Mixed LR and HR"
    })) %>%
  select(-types_list) %>%
  arrange(desc(Pattern_frequency))

# code to compile a heatmap for follow up visits
# obtain number of HPV_cols

n_follow_up <- length(HPV_cols)

# structure of matrix setting the col and row numbers and dimension name
co_occurence_fp <- matrix(" ", nrow = n_follow_up ,ncol = n_follow_up,
                          dimnames= list(HPV_cols, HPV_cols))

# Loop through lower triangle + diagonal
for (x in 1:n_follow_up){
  for (l in 1:x) {
    col_x <- as.numeric(persistent_hpv[[HPV_cols[x]]])
    col_l <- as.numeric(persistent_hpv[[HPV_cols[l]]])
    
    co_occurence_fp[x, l] <- sum(col_x == 1 & col_l == 1)}}

# Mask upper triangle with NA (to show as white)
co_occurence_fp[upper.tri(co_occurence_fp)] <- NA

# Convert to long format
heatmap_data_fp <- as.data.frame(as.table(co_occurence_fp), stringsAsFactors = FALSE)
names(heatmap_data_fp) <- c("row", "col", "value")

# convert values to numeric 
heatmap_data_fp$value <- as.numeric(heatmap_data_fp$value)

heatmap_data_fp <- heatmap_data_fp %>%
  mutate(
    row = factor(row, levels = rev(HPV_cols)),  
    col = factor(col, levels = HPV_cols))  

# PLOT --heatmap

ggplot(heatmap_data_fp, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(value), "", value)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "purple", na.value = "white") +
  coord_equal() +
  labs(x = "HPV type", y = "HPV type", fill = "Co-persistent count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("co-occurence-persist.png", width = 12, height = 6, units = "in", bg = "white")
dev.off()

# end of Objective 2 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Supplementary publication table accompanying Objective 2
# Prepare data frame from patterns dataframe by renaming column names and listing risk 
# categories using factor 

# ensure correct order of risk categories (starting with HR)
table4 <- patterns %>%
  mutate(
    risk_category_persistent = factor(
      risk_category_persistent,
      levels = c("HR only", "LR only", "Mixed LR and HR"))) %>%
  arrange(risk_category_persistent, desc(Pattern_frequency)) %>%
  rename(
    `Risk category` = risk_category_persistent,
    `Co-persistent HPV pattern` = persistent_pattern,
    `Frequency (n)` = Pattern_frequency,
    `Prevalence (%)` = percent_persistent) %>%
  select(`Risk category`, `Co-persistent HPV pattern`, `Frequency (n)`, `Prevalence (%)`)

# publication table 

table4 %>%select(-`Risk category`)%>%
  kbl(  
  caption = "Table 4: Prevalence and patterns of co-persistent HPV infection across three consecutive visits.",
  align = c("l", "r", "r"),
  booktabs = TRUE) %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12) %>%
  pack_rows("High-Risk only", 1, sum(table4$`Risk category` == "HR only"), indent = FALSE) %>%
  pack_rows("Low-Risk only",
            sum(table4$`Risk category` == "HR only") + 1,
            sum(table4$`Risk category` %in% c("HR only", "LR only")),
            indent = FALSE) %>%
  pack_rows("Mixed LR and HR",
            sum(table4$`Risk category` %in% c("HR only", "LR only")) + 1,
            nrow(table4),
            indent = FALSE) %>%
  footnote(
    general = "Total number of individuals with co-persisting infections in dataset: 64.",
    threeparttable = TRUE)

#save html table output
save_kable(table4, file = "Patterns of co-persisting infections(Obj.2).html", self_contained = TRUE)




