## Statistical analysis for PDMS bulk permeation experiment 
##
## data is provided in .csv form, 
## created by "Data Cleaning and Transformation - Bulk Diffusion.R"
##
## 31/03/2024
## Evan Kear

# Set working directory
setwd("###")

# Load packages
library(tidyverse)
library(multcomp)
library(multcompView)


# Data
data_EK <- read.csv("Performance_summary_stats_EK.csv") %>% 
  mutate(Source = "EK")
data_MB <- read.csv("Performance_summary_stats_MB.csv") %>% 
  mutate(Source = "MB")

data_time_EK <- read.csv("Performance_summary_TIME_stats_EK.csv") %>% 
  mutate(Source = "EK")
data_time_MB <- read.csv("Performance_summary_TIME_stats_MB.csv") %>% 
  mutate(Source = "MB")

data <- bind_rows(data_EK, data_MB) %>% 
  filter(Concentration == 10 | Concentration == 0)
data_time <- bind_rows(data_time_EK, data_time_MB) %>% 
  filter(Concentration == 10 | Concentration == 0)

data_300um <- bind_rows(data_EK, data_MB) %>% 
  filter(str_detect(Thickness_um, "300") | str_detect(Thickness_um, "200"))

# PDMS data
data_EK_PDMS <- data %>% 
  filter(Treatment == "PDMS" & Source == "EK") %>% 
  mutate(Thickness_um = str_replace(Thickness_um, "200 µm", "300 µm"))
data_MB_PDMS <- data %>% 
  filter(Treatment == "PDMS" & Source == "MB")


data_time_EK_PDMS <- data_time_EK %>% 
  filter(Treatment == "PDMS" & Source == "EK")
data_time_MB_PDMS <- data_time_MB %>% 
  filter(Treatment == "PDMS" & Source == "MB")


PDMS_data <- bind_rows(data_MB_PDMS, data_EK_PDMS)


# Removing duplicate PDMS and control, control from EK better, PDMS from MB used for nomalization

data <- data %>% 
  filter(!(Treatment == "PDMS" & Source == "EK")) %>% 
  filter(!(Treatment == "Control" & Source == "MB"))

data_time <- data_time %>% 
  filter(!(Treatment == "PDMS" & Source == "EK")) %>% 
  filter(!(Treatment == "Control" & Source == "MB"))

data$Treatment <- as.factor(data$Treatment)
data$Thickness_um <- as.factor(data$Thickness_um)


## Normality tests ---------------------------------------------------------
# Filter out rows where treatment or thickness is 'control'
data_filtered <- data %>%
  filter(!(Name == 'Control-0' | Thickness_um == 'control'))

# Get the unique combinations of treatment and thickness
combinations <- unique(data_filtered[, c('Name', 'Thickness_um')])


# Perform Shapiro-Wilk test for each group
shapiro_results <- data_filtered %>%
  group_by(Name, Thickness_um) %>%
  summarise(
    n = n(), # Count number of observations in the group
    # Perform test only if n >= 3
    shapiro_w_statistic = if (n() >= 3) shapiro.test(Performance_AverageFlux)$statistic else NA_real_,
    shapiro_p_value = if (n() >= 3) shapiro.test(Performance_AverageFlux)$p.value else NA_real_,
    .groups = 'drop' # Drop grouping structure after summarising
  ) %>%
  # Add an interpretation column (optional)
  mutate(
    interpretation = case_when(
      n < 3 ~ "Too few samples (n < 3)",
      shapiro_p_value > 0.05 ~ "Data likely normal (p > 0.05)",
      shapiro_p_value <= 0.05 ~ "Data likely not normal (p <= 0.05)",
      TRUE ~ "Test not performed" # Fallback for NA p-values if n>=3 (shouldn't happen)
    )
  )

# Print the results
print(shapiro_results, n=Inf) # n=Inf ensures all rows are printed

## output_filename <- "shapiro_results.txt"

# # Save the data frame to a tab-separated text file
# write.table(
#   x = shapiro_results,       # The data frame you want to save
#   file = output_filename,    # The name (and optionally path) of the output file
#   sep = "\t",                # Specify the separator (tab character)
#   row.names = FALSE,         # Don't include row numbers in the output file
#   col.names = TRUE,          # Include the column headers
#   quote = FALSE
# )


## >>>> Figure 1 C <<<< ----------------------------------
## Excipient Thickness Dunnetts tests ------------------------------------------------

data$Treatment <- as.factor(data$Treatment)
levels(data$Treatment)
data$Treatment <- relevel(data$Treatment, ref = "PDMS") # Setting PDMS as the reference for statistical analysis

data_100 <- data %>% 
  filter(Thickness_um == "100 µm")

# Fit a linear model (PDMS as control)
lm_100 <- lm(Performance_AverageFlux ~ Treatment, data = data_100)

# Run Dunnett's test (PDMS as reference level)
dunnett_test_100 <- glht(lm_100, linfct = mcp(Treatment = "Dunnett"))

# Print results
summary(dunnett_test_100)



data_300 <- data %>% 
  filter(Thickness_um == "300 µm" | Thickness_um == "200 µm")

# Fit a linear model (PDMS as control)
lm_300 <- lm(Performance_AverageFlux ~ Treatment, data = data_300)

# Run Dunnett's test (PDMS as reference level)
dunnett_test_300 <- glht(lm_300, linfct = mcp(Treatment = "Dunnett"))

# Print results
summary(dunnett_test_300)


data_500 <- data %>% 
  filter(Thickness_um == "500 µm")

# Fit a linear model (PDMS as control)
lm_500 <- lm(Performance_AverageFlux ~ Treatment, data = data_500)

# Run Dunnett's test (PDMS as reference level)
dunnett_test_500 <- glht(lm_500, linfct = mcp(Treatment = "Dunnett"))

# Print results
summary(dunnett_test_500)


## Excipient Thickness Dunnetts Loop ------------------------------------------------------------

# Ensure Treatment is a factor and set the reference level
data$Treatment <- as.factor(data$Treatment)
if ("PDMS" %in% levels(data$Treatment)) {
  data$Treatment <- relevel(data$Treatment, ref = "PDMS")
} else {
  warning("Reference level 'PDMS' not found in Treatment. Using default reference.")
}
levels(data$Treatment) # Check levels and reference

# --- Analysis for 100 µm ---
data_100 <- data %>%
  filter(Thickness_um == "100 µm")

# Check if data exists for this group before proceeding
if(nrow(data_100) > 0 && length(unique(data_100$Treatment)) > 1) {
  lm_100 <- lm(Performance_AverageFlux ~ Treatment, data = data_100)
  dunnett_test_100 <- glht(lm_100, linfct = mcp(Treatment = "Dunnett"))
  summary_100 <- summary(dunnett_test_100)
  print("--- Dunnett Test Results for 100 µm ---")
  print(summary_100) # Still print to console if desired
  
  # Extract results into a data frame
  results_100 <- data.frame(
    Thickness_Group = "100 µm",
    Comparison = names(summary_100$test$coefficients),
    Estimate = summary_100$test$coefficients,
    Std.Error = summary_100$test$sigma,
    t.value = summary_100$test$tstat,
    Adj.P.Value = summary_100$test$pvalues,
    stringsAsFactors = FALSE,
    row.names = NULL # Prevent row names based on comparison names
  )
  # Add significance stars (optional, based on common conventions)
  results_100$Signif <- symnum(results_100$Adj.P.Value, corr = FALSE, na = FALSE,
                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                               symbols = c("***", "**", "*", ".", " "))
  
} else {
  warning("Insufficient data or factor levels for 100 µm group. Skipping analysis.")
  results_100 <- NULL # Ensure variable exists but is empty
}


# --- Analysis for 200/300 µm ---
data_300 <- data %>%
  filter(Thickness_um == "300 µm" | Thickness_um == "200 µm") # Combined group

# Check if data exists for this group before proceeding
if(nrow(data_300) > 0 && length(unique(data_300$Treatment)) > 1) {
  lm_300 <- lm(Performance_AverageFlux ~ Treatment, data = data_300)
  dunnett_test_300 <- glht(lm_300, linfct = mcp(Treatment = "Dunnett"))
  summary_300 <- summary(dunnett_test_300)
  print("--- Dunnett Test Results for 200/300 µm ---")
  print(summary_300)
  
  # Extract results into a data frame
  results_300 <- data.frame(
    Thickness_Group = "200/300 µm", # Label for the combined group
    Comparison = names(summary_300$test$coefficients),
    Estimate = summary_300$test$coefficients,
    Std.Error = summary_300$test$sigma,
    t.value = summary_300$test$tstat,
    Adj.P.Value = summary_300$test$pvalues,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # Add significance stars
  results_300$Signif <- symnum(results_300$Adj.P.Value, corr = FALSE, na = FALSE,
                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                               symbols = c("***", "**", "*", ".", " "))
  
} else {
  warning("Insufficient data or factor levels for 200/300 µm group. Skipping analysis.")
  results_300 <- NULL
}


# --- Analysis for 500 µm ---
data_500 <- data %>%
  filter(Thickness_um == "500 µm")

# Check if data exists for this group before proceeding
if(nrow(data_500) > 0 && length(unique(data_500$Treatment)) > 1) {
  lm_500 <- lm(Performance_AverageFlux ~ Treatment, data = data_500)
  dunnett_test_500 <- glht(lm_500, linfct = mcp(Treatment = "Dunnett"))
  summary_500 <- summary(dunnett_test_500)
  print("--- Dunnett Test Results for 500 µm ---")
  print(summary_500)
  
  # Extract results into a data frame
  results_500 <- data.frame(
    Thickness_Group = "500 µm",
    Comparison = names(summary_500$test$coefficients),
    Estimate = summary_500$test$coefficients,
    Std.Error = summary_500$test$sigma,
    t.value = summary_500$test$tstat,
    Adj.P.Value = summary_500$test$pvalues,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # Add significance stars
  results_500$Signif <- symnum(results_500$Adj.P.Value, corr = FALSE, na = FALSE,
                               cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                               symbols = c("***", "**", "*", ".", " "))
  
} else {
  warning("Insufficient data or factor levels for 500 µm group. Skipping analysis.")
  results_500 <- NULL
}

# --- Combine all results ---
# # Use bind_rows which handles NULL inputs gracefully
# all_dunnett_results <- bind_rows(results_100, results_300, results_500)
# 
# # --- Save the combined results to a text file ---
# if (nrow(all_dunnett_results) > 0) {
#   output_filename <- "dunnett_test_summary.txt" # Choose your filename
#   
#   write.table(
#     x = all_dunnett_results,
#     file = output_filename,
#     sep = "\t",          # Tab-separated
#     row.names = FALSE,   # Don't write row numbers
#     col.names = TRUE,    # Write column headers
#     quote = FALSE        # Optional: write without quotes for cleaner look
#     # (Use TRUE if comparisons/groups might contain tabs)
#   )
#   
#   print(paste("Combined Dunnett's test results saved to:", output_filename))
#   print("--- Combined Results Data Frame ---")
#   print(all_dunnett_results) # Also print the combined table
#   
# } else {
#   print("No results were generated to save.")
# }


## Excipient Thickness Tukey's HSD ----------------------------------------------------

# Run Tukeys test
tukey_test_100 <- glht(lm_100, linfct = mcp(Treatment = "Tukey"))

# View results
summary(tukey_test_100)

# Get confidence intervals for pairwise comparisons
confint(tukey_test_100)

CLD_100 <- cld(tukey_test_100)


# Run Tukeys test
tukey_test_300 <- glht(lm_300, linfct = mcp(Treatment = "Tukey"))

# View results
summary(tukey_test_300)

# Get confidence intervals for pairwise comparisons
confint(tukey_test_300)

CLD_300 <- cld(tukey_test_300)

# Run Tukeys test
tukey_test_500 <- glht(lm_500, linfct = mcp(Treatment = "Tukey"))

# View results
summary(tukey_test_500)

# Get confidence intervals for pairwise comparisons
confint(tukey_test_500)

CLD_500 <- cld(tukey_test_500)

# # Write to text file
# sink("CLD_results.txt")
# cat("CLD for 100 μm:\n")
# print(CLD_100)
# cat("\nCLD for 300 μm:\n")
# print(CLD_300)
# cat("\nCLD for 500 μm:\n")
# print(CLD_500)
# sink()

## Excipient Thickness Tukey's HSD Loop -----------------------------------------------------------
# Function to perform Tukey's test and save results to a file
run_tukey_and_save <- function(lm_model, model_name, output_file) {
  # Run Tukey's test
  tukey_test <- glht(lm_model, linfct = mcp(Treatment = "Tukey"))
  
  # Capture summary results
  summary_results <- capture.output(summary(tukey_test))
  
  # Capture confidence intervals
  confint_results <- capture.output(confint(tukey_test))
  
  # Write results to the output file
  cat(paste("############################################################\n"), file = output_file, append = TRUE)
  cat(paste("Results for:", model_name, "\n"), file = output_file, append = TRUE)
  cat(paste("############################################################\n"), file = output_file, append = TRUE)
  
  cat("\nSummary of Tukey's Post-Hoc Test:\n", file = output_file, append = TRUE)
  cat(summary_results, sep = "\n", file = output_file, append = TRUE)
  
  cat("\nConfidence Intervals for Pairwise Comparisons:\n", file = output_file, append = TRUE)
  cat(confint_results, sep = "\n", file = output_file, append = TRUE)
  
  cat("\n\n", file = output_file, append = TRUE) # Add extra space for readability
}

# Define the output file name
output_file <- "tukey_test_results.txt"

# Create an empty output file if it doesn't exist (or overwrite if it does)
file.create(output_file)

# Run Tukey's test for lm_100 and save results
run_tukey_and_save(lm_100, "lm_100", output_file)

# Run Tukey's test for lm_300 and save results
run_tukey_and_save(lm_300, "lm_300", output_file)

# Run Tukey's test for lm_500 and save results
run_tukey_and_save(lm_500, "lm_500", output_file)




## >>>> Figure 2 B <<<< ------------------------------
## Excipient Concentration Tukey's HSD Loop  -------------------------------
data_300um$Concentration <- as.factor(data_300um$Concentration)
levels(data_300um$Concentration)
data_300um$Concentration <- relevel(data_300um$Concentration, ref = "0") # Setting PDMS as the reference for statistical analysis

# Carbopol Data
data_carbopol <- data_300um %>% 
  filter(Treatment == "Carbopol" | Treatment == "PDMS")

# Fit a linear model (PDMS as control)
lm_carbopol <- lm(Performance_AverageFlux ~ Concentration, data = data_carbopol)

# Pemulen Data
data_Pemulen <- data_300um %>% 
  filter(Treatment == "Pemulen" | Treatment == "PDMS")

# Fit a linear model (PDMS as control)
lm_Pemulen <- lm(Performance_AverageFlux ~ Concentration, data = data_Pemulen)

# PVP Data
data_PVP <- data_300um %>% 
  filter(Treatment == "PVP" | Treatment == "PDMS")

# Fit a linear model (PDMS as control)
lm_PVP <- lm(Performance_AverageFlux ~ Concentration, data = data_PVP)

# Cellulose Data
data_Cellulose <- data_300um %>% 
  filter(Treatment == "Cellulose" | Treatment == "PDMS")

# Fit a linear model (PDMS as control)
lm_Cellulose <- lm(Performance_AverageFlux ~ Concentration, data = data_Cellulose)

# CNC Data
data_CNC <- data_300um %>% 
  filter(Treatment == "CNC" | Treatment == "PDMS")

# Fit a linear model (PDMS as control)
lm_CNC <- lm(Performance_AverageFlux ~ Concentration, data = data_CNC)

# Function to perform Tukey's test and save results to a file
run_tukey_and_save_CONC <- function(lm_model, model_name, output_file_CONC) {
  # Run Tukey's test
  tukey_test <- glht(lm_model, linfct = mcp(Concentration = "Tukey"))
  
  # Capture summary results
  summary_results <- capture.output(summary(tukey_test))
  
  # Capture confidence intervals
  confint_results <- capture.output(confint(tukey_test))
  
  # CLD Calculation
  CLD_results <- cld(tukey_test)
  cld_text <- capture.output(print(CLD_results))
  
  # Write results to the output file
  cat(paste("############################################################\n"), file = output_file_CONC, append = TRUE)
  cat(paste("Results for:", model_name, "\n"), file = output_file_CONC, append = TRUE)
  cat(paste("############################################################\n"), file = output_file_CONC, append = TRUE)
  
  cat("\nSummary of Tukey's Post-Hoc Test:\n", file = output_file_CONC, append = TRUE)
  cat(summary_results, sep = "\n", file = output_file_CONC, append = TRUE)
  
  cat("\nCLD Results:\n", file = output_file_CONC, append = TRUE)
  cat(cld_text, sep = "\n", file = output_file_CONC, append = TRUE)
  
  cat("\nConfidence Intervals for Pairwise Comparisons:\n", file = output_file_CONC, append = TRUE)
  cat(confint_results, sep = "\n", file = output_file_CONC, append = TRUE)
  
  cat("\n\n", file = output_file_CONC, append = TRUE) # Add extra space for readability
}



# Define the output file name
output_file_CONC <- "tukey_test_results_CONC.txt"

# Create an empty output file if it doesn't exist (or overwrite if it does)
file.create(output_file_CONC)

# Run Tukey's test for lm_carbopol and save results
run_tukey_and_save_CONC(lm_carbopol, "lm_carbopol", output_file_CONC)

# Run Tukey's test for lm_Pemulen and save results
run_tukey_and_save_CONC(lm_Pemulen, "lm_Pemulen", output_file_CONC)

# Run Tukey's test for lm_PVP and save results
run_tukey_and_save_CONC(lm_PVP, "lm_PVP", output_file_CONC)

# Run Tukey's test for lm_Cellulose and save results
run_tukey_and_save_CONC(lm_Cellulose, "lm_Cellulose", output_file_CONC)

# Run Tukey's test for lm_CNC and save results
run_tukey_and_save_CONC(lm_CNC, "lm_CNC", output_file_CONC)




## >>>> Figure 1 B <<<< ----------------------------------------------------
# Creating subsets for all excipients, with only the 10% concentrations

## Pemulen 
data_Pemulen_10 <- data %>% 
  filter(Treatment == "Pemulen")

# Fit a linear model
lm_Pemulen_10 <- lm(Performance_AverageFlux ~ Thickness_um, data = data_Pemulen_10)

## Carbopol 
data_Carbopol_10 <- data %>% 
  filter(Treatment == "Carbopol")

# Fit a linear model
lm_Carbopol_10 <- lm(Performance_AverageFlux ~ Thickness_um, data = data_Carbopol_10)

## PVP 
data_PVP_10 <- data %>% 
  filter(Treatment == "PVP")

# Fit a linear model
lm_PVP_10 <- lm(Performance_AverageFlux ~ Thickness_um, data = data_PVP_10)

## Cellulose 
data_Cellulose_10 <- data %>% 
  filter(Treatment == "Cellulose")

# Fit a linear model
lm_Cellulose_10 <- lm(Performance_AverageFlux ~ Thickness_um, data = data_Cellulose_10)

## CNC 
data_CNC_10 <- data %>% 
  filter(Treatment == "CNC")

# Fit a linear model
lm_CNC_10 <- lm(Performance_AverageFlux ~ Thickness_um, data = data_CNC_10)

# Function to perform Tukey's test and save results to a file
run_tukey_and_save_thickness <- function(lm_model, model_name, output_file) {
  # Run Tukey's test
  tukey_test <- glht(lm_model, linfct = mcp(Thickness_um = "Tukey"))
  
  # Capture summary results
  summary_results <- capture.output(summary(tukey_test))
  
  # CLD Calculation
  CLD_results <- cld(tukey_test)
  cld_text <- capture.output(print(CLD_results))
  
  # Capture confidence intervals
  confint_results <- capture.output(confint(tukey_test))
  
  # Write results to the output file
  cat(paste("############################################################\n"), file = output_file, append = TRUE)
  cat(paste("Results for:", model_name, "\n"), file = output_file, append = TRUE)
  cat(paste("############################################################\n"), file = output_file, append = TRUE)
  
  cat("\nSummary of Tukey's Post-Hoc Test:\n", file = output_file, append = TRUE)
  cat(summary_results, sep = "\n", file = output_file, append = TRUE)
  
  
  cat("\nCLD Results:\n", file = output_file, append = TRUE)
  cat(cld_text, sep = "\n", file = output_file, append = TRUE)
  
  cat("\nConfidence Intervals for Pairwise Comparisons:\n", file = output_file, append = TRUE)
  cat(confint_results, sep = "\n", file = output_file, append = TRUE)
  
  cat("\n\n", file = output_file, append = TRUE) # Add extra space for readability
}

# Define the output file name
output_file_thickness <- "tukey_test_results_thickness.txt"

# Create an empty output file if it doesn't exist (or overwrite if it does)
file.create(output_file_thickness)

# Run Tukey's test for lm_Pemulen and save results
run_tukey_and_save_thickness(lm_Pemulen_10, "lm_Pemulen_10", output_file_thickness)

# Run Tukey's test for lm_Carbopol and save results
run_tukey_and_save_thickness(lm_Carbopol_10, "lm_Carbopol_10", output_file_thickness)

# Run Tukey's test for lm_PVP and save results
run_tukey_and_save_thickness(lm_PVP_10, "lm_PVP_10", output_file_thickness)

# Run Tukey's test for lm_Cellulose and save results
run_tukey_and_save_thickness(lm_Cellulose_10, "lm_Cellulose_10", output_file_thickness)

# Run Tukey's test for lm_CNC and save results
run_tukey_and_save_thickness(lm_CNC_10, "lm_CNC_10", output_file_thickness)


## Creating PDMS dataset
PDMS_data <- data %>% 
  filter(Treatment == "PDMS" & Source == "MB")

# Fit a linear model
lm_PDMS_10 <- lm(AverageFlux ~ Thickness_um, data = PDMS_data)

# Run Tukey's test for lm_CNC and save results
run_tukey_and_save_thickness(lm_PDMS_10, "lm_PDMS_10", output_file_thickness)
