#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Objective 1 :To determine the prevalence of multi-type HPV infections and 
#              identify which types tend to co-occur   

# Definition:
#                           
# 1. Multi-type HPV infection is infection with two or more HPV genotypes 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Libraries to load for OBJECTIVE 1
#-------------------------------------------------------------------------------
# Set up of libraries needed
# To download type install.packages("") name of library as below in "quotes"

library(tidyverse)   # collection of tools for data manipulation and visualization
library(scales)      # for axis formatting, color scales 
library(ggnewscale)  # for multiple fill/color scales in ggplot
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
# thus HPV data derived dataset
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
# Prepare data frame for baseline visits 
 
all_baseline_visits <- df_04 %>%
  filter(visit_num == 1) %>%
  select(
    study_id,
    all_of(HPV_cols),
    HPV_group,
    HPV_sum)

# Objective breakdown - Part 1: Assessing the prevalence of multi-type infections 
# pie chart of prevalence assessing distribution of Single vs Multiple vs No infection 

pie_chart_data <- all_baseline_visits %>%
  mutate(
    HPV_group = factor(HPV_group,
      levels = c("No infection", "Single infection", "Multiple infection"))) %>%
  count(HPV_group, name = "tot_individuals") %>%
  mutate(prop = round(tot_individuals / sum(tot_individuals) * 100, 1),
    legend_label = paste0(HPV_group, " ", tot_individuals, " (", prop, "%)"))

#PLOT --pie chart 

ggplot(pie_chart_data,
       aes(x = "", y = prop, fill = HPV_group)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  #labs(title = "Prevalence of HPV Infection Categories") +
  scale_fill_brewer(
    palette = "Purples",
    direction = 1,     
    name = "HPV Categories",
    labels = pie_chart_data$legend_label) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 13))


# save plot and clear plot environment
ggsave("HPV_pie_chart.png", width = 12, height = 6, units = "in" ,bg = "white")
dev.off()

# Function to create summary for stacked bar chart -clearer picture into multiple 
# infections.

hpv_summary <- function(hpv_df) {
  # total population size
  total_n <- nrow(hpv_df)
  
  results <- lapply(HPV_cols, function(hpv) {
    # subset rows where this HPV type is present
    hpv_cases <- hpv_df %>%
      filter(.data[[hpv]] == 1) %>%
      mutate(category = case_when(
        HPV_sum == 1 ~ "Single HPV genotype",
        HPV_sum == 2 ~ "Double HPV genotypes",
        HPV_sum == 3 ~ "Triple HPV genotypes",
        HPV_sum >= 4 ~ "Quad and more HPV genotypes"))
    
    # count categories and calculate population prevalence
    counts <- hpv_cases %>%
      count(category) %>%
      mutate(
        HPV_type = hpv,
        PopulationPrev = n / total_n * 100)
    
    return(counts)
  })
  
  # combine results for all HPV types
  results_df <- bind_rows(results) %>%
    rename(Count = n) %>%
    mutate(
      category = factor(
        category,
        levels = c(
          "Single HPV genotype", 
          "Double HPV genotypes", 
          "Triple HPV genotypes", 
          "Quad and more HPV genotypes"),ordered = TRUE))
  
  return(results_df)
}
# dataframe stores study_id and hpv_genotypes ONLY
# apply function to this data frame

all_hpv_genotypes <- all_baseline_visits %>%
  select(study_id,
         all_of(HPV_cols),
         HPV_sum)

#Apply function to above dataframe
summary_table <- hpv_summary(all_hpv_genotypes)
#write.csv(summary_table, "summary_table.csv", row.names = FALSE)

# Plot: Stacked bar chart of HPV genotypes by risk groups (HR/LR)
# Reorder HPV_type by total prevalence (sum across all categories), highest first

# Add risk group and calculate order

stacked_graph_data <- summary_table %>%
  group_by(HPV_type) %>%
  mutate(TotalPrev = sum(PopulationPrev)) %>%
  ungroup() %>%
  mutate(
    HPV_type = fct_reorder(HPV_type, TotalPrev, .desc = TRUE),
    risk_group = case_when(
      HPV_type %in% hr_hpv_types ~ "High-Risk",
      HPV_type %in% lr_hpv_types  ~ "Low-Risk")) 


# Data frame with total prevalence per HPV type (storing for labels)

totals <- stacked_graph_data %>%
  group_by(HPV_type, risk_group) %>%
  summarise(TotalPrev = round(sum(PopulationPrev), 1), .groups = "drop")

# PLOT --stacked bar graph
# Warnings occur because certain HPV genotypes lack double, triple, or quadruple infection categories.
# These are expected and harmless, so warnings are suppressed for cleaner output.

stacked_bar <- ggplot(stacked_graph_data, aes(x = HPV_type, y = PopulationPrev)) +
  # High-Risk bars
  geom_bar(data = filter(stacked_graph_data, risk_group == "High-Risk", na.rm = TRUE),
           aes(fill = category), stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Reds", name = "High-Risk Categories") +
  new_scale_fill() +   # reset the fill scale for Low-Risk
  
  # Low-Risk bars
  geom_bar(data = filter(stacked_graph_data, risk_group == "Low-Risk", na.rm = TRUE),
           aes(fill = category), stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Blues", name = "Low-Risk Categories") +
  
  # adding labels
  geom_text(data = totals,
            aes(x = HPV_type, y = TotalPrev, label = paste0(sprintf("%.1f", TotalPrev), "%")),
            vjust = -0.5, size = 2.5, fontface= "bold") +
  
  # Facet by risk group
  facet_wrap(~ risk_group, scales = "free_x", nrow = 1) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  labs(#title = "HPV Genotypes by Infection Category",
       x = "HPV Genotype", y = "Prevalence (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

suppressWarnings(print(stacked_bar))

# save plot and clear plot environment
suppressWarnings(
  ggsave("HPV_stacked_bar.png", width = 12, height = 6, units = "in", bg = "white"))
dev.off()

# Part 2 - co-occurrence 

# Compile a counting matrix where we see how many individuals in study have the 
# same pair of HPV. 

# vector length

n_baseline <- length(HPV_cols)

# structure of matrix setting the col and row numbers and dimension name
co_occurence <- matrix(" ", nrow = n_baseline ,ncol = n_baseline,
                       dimnames= list(HPV_cols, HPV_cols))

# Loop through lower triangle + diagonal
for (k in 1:n_baseline){
  for (j in 1:k) {   # j starts at k fills lower triangle
    col_k <- as.numeric(all_baseline_visits[[ HPV_cols[k] ]])
    col_j <- as.numeric(all_baseline_visits[[ HPV_cols[j] ]])
    
    co_occurence[k, j] <- sum(col_k == 1 & col_j == 1)}}

# Mask upper triangle with NA (to show as white)
co_occurence[upper.tri(co_occurence)] <- NA

# Convert to long format
heatmap_data <- as.data.frame(as.table(co_occurence), stringsAsFactors = FALSE)
names(heatmap_data) <- c("row", "col", "value")

# convert values to numeric 
heatmap_data$value <- as.numeric(heatmap_data$value)

heatmap_data <- heatmap_data %>%
  mutate(
    # reverse order for y-axis
    row = factor(row, levels = rev(HPV_cols)),  
    col = factor(col, levels = HPV_cols))       

# PLOT --heatmap 
ggplot(heatmap_data, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(value), "", value)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "purple", na.value = "white") +
  coord_equal() +
  labs(x = "HPV type", y = "HPV type", fill = "Dual infecion count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# save plot and clear plot environment
ggsave("co-occurence.png", width = 12, height = 6, units = "in", bg = "white")
dev.off()
#================================================================================
#================================================================================
# Exact pattern co-occurrence(identifying each individuals infection profile)
# Create a pattern column per individual and reporting the pattern using collapse

# dataframe stores baseline visits AND multiple HPV infected individuals
# case where individuals are infected with more than one HPV genotype

baseline_visits_multiple <- all_baseline_visits %>%
  filter(HPV_group  == "Multiple infection" ) %>%
  select(study_id,
         (all_of(HPV_cols)))

individual_pattern <- baseline_visits_multiple %>%
  mutate(
    pattern = apply(select(., all_of(HPV_cols)), 1, function(x) {
      pos <- HPV_cols[which(x == 1)]
      if (length(pos) == 0) "No infection"
      else paste(pos, collapse = " + ")
    }))

# Add common infection profiles
hpv_pattern <- individual_pattern %>%
  select(study_id, pattern) %>%
  count(pattern, name = "Pattern_count") %>%
  mutate(percent = round((Pattern_count / sum(Pattern_count)) * 100, 2))%>%
  # Classify risk category for each pattern
  mutate(
    types_list = str_split(pattern, " \\+ ") %>% map(~ str_trim(.x)), 
    risk_category = map_chr(types_list, function(types) {
      if (all(types %in% hr_hpv_types)) {
        "HR only"
      } else if (all(types %in% lr_hpv_types)) {
        "LR only"
      } else {
        "Mixed LR and HR"
      }
    })
  ) %>%
  select(-types_list)%>%
  # arrange list in descending order
  arrange(desc(Pattern_count))

# end of objective 1 
#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------
# # Publication style tables 
# supplementary table for Objective 1
# For each HPV genotype, counts show the number (and percentage) of samples 
# where that genotype occurred by itself (single) or concurrently 

# obtain unique individual number in data frame
n_total <- nrow(all_baseline_visits)

# Recode categories into single vs multiple ONLY (same as above but limiting to single VS multiple)
# mutate column that classifies single and multiple infection 

summary_new <- summary_table %>%
  mutate(
    infection_group = ifelse(category == "Single HPV genotype", "Single_infection", "Multiple_infection"))

# Summarize counts by HPV type and infection group
# and pivot data frame to wide format

summary_counts <- summary_new %>%
  group_by(HPV_type, infection_group) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    PopulationPrev = sum(PopulationPrev, na.rm = TRUE),
    .groups = "drop")%>%
  pivot_wider(
      names_from = infection_group,
      values_from = c(Count, PopulationPrev),
      names_glue = "{.value}_{infection_group}",
      values_fill = 0)

wide_format <- summary_counts %>%
  mutate(
    Total_Positive = Count_Single_infection + Count_Multiple_infection,
    Positive_Percent = round((Total_Positive / n_total) * 100, 2))%>%
# Compute 95% CI for each HPV type 
  rowwise() %>%
  mutate(
    ci = list(broom::tidy(prop.test(Total_Positive, n_total)) %>%
                select(conf.low, conf.high) %>%
                mutate(across(everything(), ~ round(.x * 100, 2)))),
    CI_95 = paste0(ci$conf.low, "–", ci$conf.high)) %>%
  select(-ci) %>%
  ungroup()


# prepare data frame for table by renaming col names and ensuring ordering of 
# HPV genotypes from high risk to low risk

table2 <- wide_format %>%
  transmute(
    `HPV Genotype` = HPV_type,
    `Positive n(%)` = paste0(Total_Positive, " (", Positive_Percent, ")"),
    `95% CI for positive` = CI_95,
    `Single-type infection n(%)` = paste0(Count_Single_infection, " (", round((Count_Single_infection / n_total) * 100, 2), ")"),
    `Multiple-tpe infection n(%)` = paste0(Count_Multiple_infection, " (", round((Count_Multiple_infection / n_total) * 100, 2), ")")) %>%
  mutate(
    # classify for table ordering 
    risk_order = case_when(
      `HPV Genotype` %in% hr_hpv_types ~ 1,
      `HPV Genotype` %in% lr_hpv_types  ~ 2,),
    numeric_type = as.numeric(gsub("HPV\\d+", "", `HPV Genotype`))) %>%
  arrange(risk_order, numeric_type) %>%
  select(-risk_order, -numeric_type)


# publication-style table 
table2 %>%
  kbl(caption = "Table 2: Distribution of single vs multiple HPV genotype infections.",
    align = c("l", "r", "r", "r", "r"),
    booktabs = TRUE,
    escape = FALSE) %>%
  kable_styling(
    full_width = FALSE,position = "center",font_size = 12) %>%
  pack_rows(
    "High-Risk only",1,sum(table2$`HPV Genotype` %in% hr_hpv_types),indent = FALSE) %>%
  pack_rows(
    "Low-Risk only",sum(table2$`HPV Genotype` %in% hr_hpv_types) + 1,nrow(table2),
    indent = FALSE)

# save html output table 

save_kable(table2, file = "Single Vs Multiple infections(Obj.1).html",self_contained = TRUE)
#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------
# Publication style table 
# Supplementary table 3

# Full detailed table of each individuals infection profile  
# ensure correct order of risk categories
table3 <- hpv_pattern %>%
  mutate(
    risk_category = factor(
      risk_category,
      levels = c("HR only", "LR only", "Mixed LR and HR"))) %>%
  #arrange(risk_category, desc(Pattern_count)) %>%
  rename(
    `Risk category` = risk_category,
    `HPV pattern` = pattern,
    `Frequency (n)` = Pattern_count,
    `Prevalence (%)` = percent) %>%
  select(`Risk category`, `HPV pattern`, `Frequency (n)`, `Prevalence (%)`)

# publication table 

table3 %>% select(-`Risk category`)%>%
  kbl(
  caption = "Table 3: HPV infection patterns and prevalence at baseline.",
  align = c("l", "r", "r"),
  booktabs = TRUE) %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12) %>%
  pack_rows(
    "High-Risk only",1,sum(table3$`Risk category` == "HR only"),
    indent = FALSE) %>%
  pack_rows(
    "Low-Risk only",
    sum(table3$`Risk category` == "HR only") + 1,
    sum(table3$`Risk category` %in% c("HR only", "LR only")),
    indent = FALSE) %>%
  pack_rows(
    "Mixed LR and HR",
    sum(table3$`Risk category` %in% c("HR only", "LR only")) + 1,nrow(table3),
    indent = FALSE)

# save html output table 
save_kable(table3, file = "HPV patterns(Obj.1)",self_contained = TRUE)

