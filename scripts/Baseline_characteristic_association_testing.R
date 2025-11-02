#-------------------------------------------------------------------------------
# Baseline characteristic hypothesis testing
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Set up of libraries needed
# To download type install.packages("") name of library as below in "quotes"

library(tidyverse) # collection of tools for data manipulation and visualization
library(broom)     # p-value formatting 
library(knitr)     # for knitting tables/reports
library(kableExtra)# for publication-ready tables
library(rlang)     # for tidy evaluation and programming utilities

# Load HPV genotyping containing dataset 

HPV_genotyping_data <- read_csv("C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts\\HPV_genotyping_complete.csv")
metadata_clean_03 <- read_csv("C:\\Users\\khomo\\Downloads\\Honours_Project_ 2025\\HPV_raw_data_and_scripts\\Sexual_behaviour_metadata_final_clean.csv")

# Keep only baseline visits (visit_num == 1) and  

all_baseline_visits <- HPV_genotyping_data %>%
  filter(visit_num == 1)

#-------------------------------------------------------------------------------
# keep only individuals who appear in both data frames via inner_join
# thus everyone in baseline_data that has HPV info and metadata

baseline_data <- metadata_clean_03 %>%
  inner_join(all_baseline_visits, by = "study_id")

# vector containing HPV columns 
HPV_cols <- grep("HPV_\\d+", names(baseline_data), value= TRUE)

# encode numeric column HPV_sum into HPV- and HPV+ Single, HPV+ Multiple
baseline_data <- baseline_data %>%
  mutate(
    infection_status = case_when(
      HPV_sum == 0 ~ "HPV-negative",
      HPV_sum == 1 ~ "Single infection",
      HPV_sum > 1  ~ "Multiple infection",
      TRUE ~ NA_character_),
    HPV_group = case_when(
      infection_status == "HPV-negative" ~ "HPV−",
      infection_status == "Single infection" ~ "HPV+ Single",
      infection_status == "Multiple infection" ~ "HPV+ Multiple"),
    HPV_group = factor(
      HPV_group,
      levels = c("HPV−", "HPV+ Single", "HPV+ Multiple")))

#-------------------------------------------------------------------------------
# Recode numeric variables into descriptive categorical (factor) labels

baseline_data <- baseline_data %>%
  mutate(
    condoms_used_30d = factor(condoms_used_30d,
                              levels = c(1, 2, 3),
                              labels = c("Never", "Sometimes", "Always")),
    penile_vaginal_sex_7d = factor(penile_vaginal_sex_7d,
                                   levels = c(0, 1),
                                   labels = c("No", "Yes")),
    oral_sex_7d = factor(oral_sex_7d,
                         levels = c(0, 1),
                         labels = c("No", "Yes")),
    anal_sex_7d = factor(anal_sex_7d,
                         levels = c(0, 1),
                         labels = c("No", "Yes")),
    sti_symptoms_now = factor(sti_symptoms_now,
                              levels = c(0, 1),
                              labels = c("No", "Yes")),
    live_with_partner_baseline = factor(live_with_partner_baseline,
                                        levels = c(0, 1),
                                        labels = c("No", "Yes")))
#-------------------------------------------------------------------------------
# Function to perform one-way ANOVA for normally distributed variables 
# and generate a summary table with mean ± SD and p-value per group

anova_summary <- function(data, group_var, vars, var_labels = NULL) {
  group_var <- ensym(group_var)
  group_name <- as_string(group_var)
  
  # Count totals for column headers
  total_n <- nrow(filter(data, !is.na(!!group_var)))
  group_counts <- table(data[[group_name]])
  group_levels <- names(group_counts)
  
  # Build results table
  results <- map_dfr(vars, function(var) {
    var_sym <- ensym(var)
    var_name <- as_string(var_sym)
    display_name <- if (!is.null(var_labels) && var_name %in% names(var_labels)) {
      var_labels[[var_name]]
    } else {
      var_name
    }
    
    # Overall mean ± SD
    overall <- data %>%
      summarise(
        Mean = mean(!!var_sym, na.rm = TRUE),
        SD   = sd(!!var_sym, na.rm = TRUE)
      )
    
    # Group-wise mean ± SD
    group_stats <- data %>%
      group_by(!!group_var) %>%
      summarise(
        Mean = mean(!!var_sym, na.rm = TRUE),
        SD   = sd(!!var_sym, na.rm = TRUE),
        n    = sum(!is.na(!!var_sym)),
        .groups = "drop"
      )
    
    # Perform one-way ANOVA
    aov_res <- aov(data[[var_name]] ~ data[[group_name]])
    p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
    
    # Build formatted row
    tibble(
      Variable = paste0(display_name, " (mean ± SD)"),
      Total = sprintf("%.2f ± %.2f", overall$Mean, overall$SD),
      !!!set_names(
        map(group_levels, function(g) {
          sprintf("%.2f ± %.2f",
                  group_stats$Mean[group_stats[[group_name]] == g],
                  group_stats$SD[group_stats[[group_name]] == g])
        }),
        group_levels
      ),
      `p-value` = formatC(p_val, format = "f", digits = 3)
    )
  })
  
  # Publication-style table
  results %>%
    kbl(
      align = "lcccc",
      booktabs = TRUE,
      escape = FALSE,
      col.names = c(
        "Variable",
        sprintf("Total (n=%d)", total_n),
        sprintf("%s (n=%d)", group_levels[1], group_counts[1]),
        sprintf("%s (n=%d)", group_levels[2], group_counts[2]),
        sprintf("%s (n=%d)", group_levels[3], group_counts[3]),
        "p-value"
      )
    ) %>%
    kable_styling(full_width = FALSE, position = "center", font_size = 12)
  return(results)
}


#-------------------------------------------------------------------------------
# Function to perform Kruskal–Wallis tests for non-normally distributed variables 
# and generate a summary table with median, IQR, and p-value per group


kw_summary <- function(data, group_var, vars, var_labels = NULL) {
  group_var <- ensym(group_var)
  group_name <- as_string(group_var)
  
  # group counts
  total_n <- nrow(filter(data, !is.na(!!group_var)))
  group_counts <- table(data[[group_name]])
  group_levels <- names(group_counts)
  
  
  results <- map_dfr(vars, function(var) {
    var_sym <- ensym(var)
    var_name <- as_string(var_sym)
    display_name <- if (!is.null(var_labels) && var_name %in% names(var_labels)) {
      var_labels[[var_name]]
    } else var_name
    
    # overall median + IQR
    overall <- data %>% summarise(
      Median = median(!!var_sym, na.rm = TRUE),
      Q1 = quantile(!!var_sym, 0.25, na.rm = TRUE),
      Q3 = quantile(!!var_sym, 0.75, na.rm = TRUE)
    )
    
    # group-wise summary
    summary_stats <- data %>%
      group_by(!!group_var) %>%
      summarise(
        Median = median(!!var_sym, na.rm = TRUE),
        Q1 = quantile(!!var_sym, 0.25, na.rm = TRUE),
        Q3 = quantile(!!var_sym, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Kruskal–Wallis test (handles >2 groups)
    kw_res <- kruskal.test(data[[var_name]] ~ data[[group_name]])
    
    # formatted row
    tibble(
      Variable = paste0(display_name, " (median [IQR])"),
      Total = sprintf("%.2f [%.2f–%.2f]", overall$Median, overall$Q1, overall$Q3),
      !!!set_names(
        map(group_levels, function(g) {
          sprintf("%.2f [%.2f–%.2f]",
                  summary_stats$Median[summary_stats[[group_name]] == g],
                  summary_stats$Q1[summary_stats[[group_name]] == g],
                  summary_stats$Q3[summary_stats[[group_name]] == g])
        }),
        group_levels
      ),
      `p-value` = formatC(kw_res$p.value, format = "f", digits = 3)
    )
  })
  
  # Publication-style table
  results %>%
    kbl(
      col.names = c(
        "**Variable**",
        sprintf("**Total (n=%d)**", total_n),
        sprintf("**%s (n=%d)**", group_levels[1], group_counts[1]),
        sprintf("**%s (n=%d)**", group_levels[2], group_counts[2]),
        sprintf("**%s (n=%d)**", group_levels[3], group_counts[3]),
        "**p-value**"
      ),
      align = "lccccc",
      booktabs = TRUE,
      escape = FALSE
    ) %>%
    kable_styling(full_width = FALSE, position = "center", font_size = 12)
  
  return(results)
}


#-------------------------------------------------------------------------------
# Function to perform Chi-square or Fisher’s exact tests for categorical variables 
# and generate a summary table with counts, percentages, and p-values per group


categorical_summary <- function(data, group_var, vars, var_labels = NULL) {
  group_var <- ensym(group_var)
  group_name <- as_string(group_var)
  
  # header counts
  total_n <- nrow(filter(data, !is.na(!!group_var)))
  group_counts <- table(data[[group_name]])
  group_levels <- names(group_counts)
  
  results <- map_dfr(vars, function(var) {
    var_sym <- ensym(var)
    var_name <- as_string(var_sym)
    display_name <- if (!is.null(var_labels) && var_name %in% names(var_labels)) {
      var_labels[[var_name]]
    } else var_name
    
    # ensure factor levels preserved for consistent ordering
    df_var <- data %>%
      mutate(!!var_sym := factor(!!var_sym, levels = unique(na.omit(data[[var_name]])))) %>%
      filter(!is.na(!!var_sym))
    
    if (nrow(df_var) == 0) return(NULL)
    
    # contingency table + test
    tbl <- table(df_var[[var_name]], df_var[[group_name]])
    if (any(dim(tbl) == 0)) return(NULL)
    
    test_used <- if (any(suppressWarnings(chisq.test(tbl))$expected < 5)) "Fisher" else "Chi-square"
    p_val <- if (test_used == "Fisher") fisher.test(tbl)$p.value else suppressWarnings(chisq.test(tbl)$p.value)
    
    # overall + group counts
    total_counts <- df_var %>%
      count(!!var_sym) %>%
      mutate(Percent = (n / sum(n)) * 100)
    
    group_counts_df <- df_var %>%
      group_by(!!group_var, !!var_sym) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(!!group_var) %>%
      mutate(Percent = (n / sum(n)) * 100)
    
    # header row
    header_row <- tibble(
      Variable = display_name,
      Total = "",
      !!!set_names(rep(list(""), length(group_levels)), group_levels),
      `p-value` = paste0(formatC(p_val, format = "f", digits = 3),
                         if (test_used == "Fisher") "*" else "")
    )
    
    # category rows
    category_rows <- map_dfr(levels(df_var[[var_name]]), function(cat) {
      total <- sprintf("%d (%.1f%%)",
                       total_counts$n[total_counts[[var_name]] == cat],
                       total_counts$Percent[total_counts[[var_name]] == cat])
      
      group_values <- map_chr(group_levels, function(g) {
        n_val <- group_counts_df$n[group_counts_df[[group_name]] == g & group_counts_df[[var_name]] == cat]
        p_val <- group_counts_df$Percent[group_counts_df[[group_name]] == g & group_counts_df[[var_name]] == cat]
        if (length(n_val) == 0) return("0 (0.0%)")
        sprintf("%d (%.1f%%)", n_val, p_val)
      })
      
      tibble_row <- tibble(
        Variable = paste0("   ", as.character(cat)),
        Total = total,
        !!!set_names(as.list(group_values), group_levels),
        `p-value` = ""
      )
      tibble_row
    })
    
    bind_rows(header_row, category_rows)
  })
  
  # Final formatted publication-style table
  results %>%
    kbl(
      col.names = c(
        "**Variable**",
        sprintf("**Total (n=%d)**", total_n),
        sprintf("**%s (n=%d)**", group_levels[1], group_counts[1]),
        sprintf("**%s (n=%d)**", group_levels[2], group_counts[2]),
        sprintf("**%s (n=%d)**", group_levels[3], group_counts[3]),
        "**p-value**"
      ),
      align = "lccccc",
      booktabs = TRUE,
      escape = FALSE
    ) %>%
    kable_styling(full_width = FALSE, position = "center", font_size = 12)
  
  return(results)
}

#===============================================================================
# Compute and label ANOVA, Kruskal–Wallis, and categorical summaries comparing 
# baseline characteristics across HPV groups


vars_anova <- c("age_at_baseline", "age_first_sex_baseline")

anova_df <- anova_summary(
  baseline_data,group_var = HPV_group,vars = vars_anova,
  var_labels = c(
    "age_at_baseline" = "Age at baseline(mean ± SD)",
    "age_first_sex_baseline" = "Age at sexual debut(mean ± SD)"))%>%
  mutate(
    Variable = ifelse(`p-value` != "" & !is.na(`p-value`),
                      paste0(Variable, "<sup>1</sup>"),
                      Variable))

vars_kw <- c(
  "number_sexual_partners","sex_episodes_7d","sex_partners_7d","sex_episodes_30d",
  "sex_partners_30d")
kw_df   <- kw_summary(baseline_data, group_var = HPV_group,
                              vars = vars_kw,
                              var_labels = c(
                                "number_sexual_partners" = "Number of sexual partners(median[IQR])",
                                "sex_episodes_7d" = "Sex events in past 7days(median[IQR]) ",
                                "sex_partners_7d" = "Number of sex partners in past 7 days(median[IQR]) ",
                                "sex_episodes_30d" = "Sex events in past 30days(median[IQR]) ",
                                "sex_partners_30d" = "Number of sex partners in past 30 days(median[IQR]) ")) %>%
  mutate(
    Variable = ifelse(`p-value` != "" & !is.na(`p-value`),
                      paste0(Variable, "<sup>2</sup>"),
                      Variable))


vars_categorical <- c(
   "condoms_used_30d","sti_symptoms_now")

categorical_df <- categorical_summary(baseline_data, group_var = HPV_group, 
                                      vars = vars_categorical,
                                      var_labels = c(
                                        "condoms_used_30d" = "Condom use in past month (n(%))",
                                        "sti_symptoms_now" = "Current STI symptoms (n(%))"))%>%
  mutate(
    Variable = ifelse(`p-value` != "" & !is.na(`p-value`),
                      paste0(Variable, "<sup>3</sup>"),
                      Variable))

# combine all dataframes with bind_rows 

combined_df <- bind_rows(anova_df ,kw_df,categorical_df)

# Define header totals 

total_n <- nrow(filter(baseline_data, !is.na(HPV_group )))
n_neg    <- sum(baseline_data$HPV_group == "HPV−")
n_single <- sum(baseline_data$HPV_group == "HPV+ Single")
n_multi  <- sum(baseline_data$HPV_group == "HPV+ Multiple")


# Produce publication table

combined_df %>%
  kbl(
    caption = "Table 1: HPV Baseline Characteristics<br><span style='font-size:11px'><sup>1</sup>ANOVA test; <sup>2</sup>Kruskal–Wallis test; <sup>3</sup>Chi-square.</span><br><span style='font-size:11px'></span>",
    align = "lcccc",
    booktabs = TRUE,
    escape = FALSE,
    col.names = c(
      "Variable",
      sprintf("Total (n=%d)", total_n),
      sprintf("HPV− (n=%d)", n_neg),
      sprintf("HPV+ Single (n=%d)", n_single),
      sprintf("HPV+ Multiple (n=%d)", n_multi),
      "p-value")) %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12)


