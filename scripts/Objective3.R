#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Objective 3 : To classify individuals into distinct longitudinal infection 
#               trajectories       

# Longitudinal trajectories:
#                           
# 1. Clearing multi-type (first visit = Multi-type, last visit = Negative)
# 2. Incident multi-type (first visit = Single, last visit = Multi-type)
# 3. Persistent multi-type (ALL 3 visits individuals are infected with Multi-type)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Libraries to load for OBJECTIVE 3

#-------------------------------------------------------------------------------
# Set up of libraries needed
# To download type install.packages("") name of library as below in "quotes"

library(tidyverse)   # collection of tools for data manipulation and visualization
library(knitr)       # for knitting tables/reports
library(kableExtra)  # for publication-ready tables 
library(ggalluvial)  # for alluvial(flow) plots

#-------------------------------------------------------------------------------
#set path to working directory

WKDIR <- "C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts"
setwd(WKDIR)

# Read in csv
# first containing individuals with 3 visits 
# second persisting infection (collapsed three visits to one entry per individual)
longitudinal_trajectories <- read_csv("complete_consecutive_visits.csv")
persistent_hpv <-read_csv("persistening_hpv_patterns.csv")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# create vectors to store HPV genotype columns 
# first line is all columns (all HPV types in data frame) 
# second line is high risk types only 
# third line is low risk types only

HPV_cols <- grep("HPV\\d+", names(longitudinal_trajectories), value= TRUE)

hr_hpv_types <- c("HPV16", "HPV18", "HPV26", "HPV31", "HPV33", "HPV35",
                  "HPV39", "HPV45", "HPV51", "HPV52", "HPV53", "HPV56",
                  "HPV58", "HPV59", "HPV66", "HPV68", "HPV73", "HPV82")

lr_hpv_types <- c("HPV06", "HPV11", "HPV40", "HPV42", "HPV43",
                  "HPV54", "HPV61", "HPV67", "HPV69",
                  "HPV70", "HPV71", "HPV72", "HPV84")


#------------------------------------------------------------------------------#
# Define states(Negative/single/multiple infection) to allocate individuals visit 


states_df <- longitudinal_trajectories %>%
  # Identify positive HPV types per individual
  mutate(
    typeset = apply(select(., all_of(HPV_cols)), 1, function(x) {
      pos <- names(x)[x == 1]  # HPV types that are positive
      if (length(pos) == 0) NA_character_ else paste(sort(pos), collapse = " + ")}),
    num_types = ifelse(is.na(typeset), 0L, str_count(typeset, "\\+") + 1L),
    State = case_when(
      num_types == 0 ~ "Negative",
      num_types == 1 ~ "Single infection",
      num_types >= 2 ~ "Multiple infection")) %>%
  select(study_id, visit_num, typeset, num_types, State)


# Endpoint-based trajectories( between first and last visit)
# Classify individuals into trajectories - Clearing,Incident and Other 

endpoints_df <- states_df %>%
  arrange(study_id, visit_num) %>%
  group_by(study_id) %>%
  summarise(
    first_state = first(State),
    last_state  = last(State),
    .groups = "drop") %>%
  mutate(Trajectory = case_when(
    first_state == "Multiple infection" & last_state == "Negative" ~ "Clearing Multiple infection",
    first_state == "Single infection"   & last_state == "Multiple infection" ~ "Incident Multiple infection",
    TRUE ~ "Other"))  # individuals NOT fitting pre defined trajectories

  

# Filter table to contain those with 2 or more (multiple infection) from persisting dataframe
# and determine risk category 

multiple_persisting <- persistent_hpv %>%
  mutate(n_pos = rowSums(select(., all_of(HPV_cols)), na.rm = TRUE)) %>%
  filter(n_pos >= 2) %>%
  select(-n_pos) %>%
# determine risk category for each pattern
  mutate(
    # Identify which HPV types persisted (where value == 1)
    pattern_persistent = apply(select(., all_of(HPV_cols)), 1, function(x) {
      pos <- names(x)[x == 1]
      if (length(pos) == 0) NA_character_ else paste(pos, collapse = " + ")}),
    
    # Split and trim pattern for classification
    types_list = str_split(pattern_persistent, " \\+ ") %>% map(~ str_trim(.x)),
    
    # Determine risk classification for each pattern
    risk_category_persistent = map_chr(types_list, function(types) {
      if (all(types %in% hr_hpv_types)) {
        "HR only"
      } else if (all(types %in% lr_hpv_types)) {
        "LR only"
      } else {
        "Mixed LR and HR"
      }
    })
  ) %>%
  select(study_id,risk_category_persistent)

# Identify true persisters among the "Other" category using the study_ids from 
# multiple_persisting data frame

endpoints_reclassified <- endpoints_df %>%
  mutate(
    Trajectory = case_when(
      Trajectory == "Incident Multiple infection" ~ "Incident Multiple infection",
      Trajectory == "Clearing Multiple infection" ~ "Clearing Multiple infection",
      Trajectory == "Other" & study_id %in% multiple_persisting$study_id ~ "Persisting Multiple infection",
      TRUE ~ "Other"))

# summarise total counts 
# Count infection trajectories
# Obtain total individuals 
total_n <- nrow(endpoints_reclassified)

summary_endpoints <- endpoints_reclassified %>%
  count(Trajectory, name = "n") %>%
  mutate(percent = round((n / total_n) * 100, 1),
         label = paste0(n, " (", percent, ")")) %>%
  select(Trajectory, label)


# Obtain total risk category - HR/LR and mixed

summary_persisting <- multiple_persisting %>%
  count(risk_category_persistent, name = "n") %>%
  mutate(
    percent = round((n / total_n) * 100, 1),
    label = paste0(n, " (", percent, ")"),
    Trajectory = "Persisting Multiple infection",
    risk_category_persistent = factor(
      risk_category_persistent,
      levels = c("HR only", "LR only", "Mixed LR and HR"))) %>%
  arrange(risk_category_persistent) %>%
  filter(!is.na(risk_category_persistent)) %>%
  # Create indented rows for the persisting risk categories
  mutate(Trajectory = paste0("   ", risk_category_persistent)) %>%
  select(Trajectory, label)


# Prepare dataframe for plot into long format and use factor to ensure a ordered level

endpoints_long <- endpoints_reclassified %>%
  pivot_longer(
    cols = c(first_state, last_state),
    names_to = "VisitStage",
    values_to = "State") %>%
  mutate(
    VisitStage = recode(VisitStage,
                        "first_state" = "First visit",
                        "last_state"  = "Last visit"),
    VisitStage = factor(VisitStage, levels = c("First visit", "Last visit")),
    State = factor(State, levels = c("Negative", "Single infection", "Multiple infection")),
    # Controls legend order
    Trajectory = factor(Trajectory,
                        levels = c("Persisting Multiple infection",
                                   "Incident Multiple infection",
                                   "Clearing Multiple infection",
                                   "Other")))

# Count percentage for creating label for plots legend 
traj_percents <- endpoints_long %>%
  distinct(study_id, Trajectory) %>%    # one row per participant
  count(Trajectory) %>%
  mutate(
    perc = round(100 * n / sum(n), 1),
    label = paste0(n, " (", perc, "%)"))    

# Add label for legends created to data frame for plot 
endpoints_long <- endpoints_long %>%
  left_join(traj_percents, by = "Trajectory") %>%
  mutate(
    Trajectory_label = paste0(Trajectory, "  ", label),
    Trajectory_label = factor(
      Trajectory_label,
      levels = paste0(traj_percents$Trajectory, "  ", traj_percents$label)))  



# PLOT --alluvial

ggplot(endpoints_long,
       aes(x = VisitStage,
           stratum = State,
           alluvium = study_id,
           fill = Trajectory_label,
           order = as.numeric(Trajectory_label))) + 
  # Ribbon width made a bit wide
  geom_alluvium(alpha = 0.9, width = 0.15) +   
  geom_stratum(width = 0.30, color = "black", fill = "grey95") +
  # Labels for each stratum (Negative, Single, Multiple)
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Purples", direction = -1, name = "Trajectory",na.translate = FALSE ) +
  # Title and theme
  labs(
    title = "HPV Infection Trajectories Across First and Last Visits",
    x = "Visit Stage",
    y = "Number of Individuals") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"))

# save plot and clear plot environment
ggsave("alluvia.png", width = 12, height = 6, units = "in", bg = "white")
dev.off()

# end of Objective 3 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Supplementary publication table accompanying Objective 3
# Prepare data frame by binding rows that have summary of trajectories and counts
 
table5 <- bind_rows(
  summary_endpoints,
  summary_persisting) %>%
  arrange(match(Trajectory, c(
    "Incident Multiple infection",
    "Clearing Multiple infection",
    "Persisting Multiple infection",
    "   HR only",
    "   LR only",
    "   Mixed LR and HR",
    "Other")))


# publication table 

table5 %>% 
  kbl(caption = "Table 5: Longitudinal infection trajectories.",
    col.names = c("Trajectory type", "Frequency(%)"),
    align = c("l", "r"),
    booktabs = TRUE,
    escape = TRUE) %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12) %>%
  row_spec(3, bold = TRUE, background = "#f0f0f0") %>%   # highlight main persisting row
  footnote(
    general = paste0(
      "‘Other’ includes participants with undefined infection trajectories. ",
      "Total number of participants included in longitudinal analysis: ", total_n, "."),
    threeparttable = TRUE)

#save html table output
save_kable(table5, file = "Longitudinal_trajectories(Obj.3).html", self_contained = TRUE)

