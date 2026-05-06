## Fructose diffusion analysis: visualization, and statistics
## Inputs: Combined CSV data with YFP, and flux columns
## Outputs: Plots, normality tests, ANOVA with post-hoc comparisons
## Author: Evan Kear, Updated: 2025-11-17


# Load libraries
library(tidyverse)
library(patchwork)
library(multcomp)
library(multcompView)

# Set working directory
wd <- "XXX"
setwd(wd)

## Data prep ---------------------------------------------------------------
Data <- read.csv("All_Data.csv")

Data_filtered <- Data %>%
  filter(AUC_YFP > 0) %>%
  mutate(Log10_Flux = as.numeric(Log10_Flux))

Data_summary <- Data_filtered %>%
  group_by(Condition, Sample.Point) %>%
  summarise(
    AVG_Log_Flux = mean(Log10_Flux, na.rm = TRUE),
    SD_Log_Flux = sd(Log10_Flux, na.rm = TRUE),
    AVG_Mass_Flux = mean(Mass_Flux_kgms, na.rm = TRUE),
    SD_Mass_Flux = sd(Mass_Flux_kgms, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

Data_main <- Data_filtered %>%
  left_join(Data_summary, by = c("Condition", "Sample.Point"))

# set factors
colour_list <- c("Control" = "grey50", "PDMS" = "#A32A3A", "CNC" = "#507E4D", "Cellulose" = "#E0AE51", 
                 "PVP" = "#F58026", "Pemulen" = "#C66BF5", "Carbopol" = "#0C34AB", "Cuticle" = "#00CD00")

Data_main$Condition <- factor(Data_main$Condition, 
                              levels = c("Control", "PDMS", "CNC", "Cellulose", 
                                         "PVP", "Pemulen", "Carbopol", "Cuticle"))

# DATA FOR T1 AND T2
t1_data <- Data_main %>% filter(Sample.Point == "T1")
t2_data <- Data_main %>% filter(Sample.Point == "T2")


## Plotting ----------------------------------------------------------------
# T1 PLOTS
Log_T1 <- t1_data %>%
  ggplot(aes(x = Condition, y = Log10_Flux, colour = Condition)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_errorbar(aes(ymin = AVG_Log_Flux - SD_Log_Flux,
                    ymax = AVG_Log_Flux + SD_Log_Flux),
                width = 0.2, colour = "black") +
  geom_errorbar(
    aes(ymin = AVG_Log_Flux, ymax = AVG_Log_Flux),
    width = 0.5, colour = "black") +
  scale_y_continuous(name = expression('Log10 Flux (kg·m'^-2*'s'^-1*')'),
                     expand = expansion(mult = c(0, 0.05)),
                     limits = c(-11, -4)) +
  scale_x_discrete(name = "Treatment") +
  scale_colour_manual(values = colour_list) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))


Mass_T1 <- t1_data %>%
  ggplot(aes(x = Condition, y = Mass_Flux_kgms, colour = Condition)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_errorbar(aes(ymin = AVG_Mass_Flux - SD_Mass_Flux,
                    ymax = AVG_Mass_Flux + SD_Mass_Flux),
                width = 0.2, colour = "black") +
  geom_errorbar(
    aes(ymin = AVG_Mass_Flux, ymax = AVG_Mass_Flux),
    width = 0.5, colour = "black") +
  scale_y_continuous(name = expression('Flux (kg·m'^-2*'s'^-1*')'),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(name = "Treatment") +
  scale_colour_manual(values = colour_list) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

# T2 PLOTS
Log_T2 <- t2_data %>%
  ggplot(aes(x = Condition, y = Log10_Flux, colour = Condition)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_errorbar(aes(ymin = AVG_Log_Flux - SD_Log_Flux,
                    ymax = AVG_Log_Flux + SD_Log_Flux),
                width = 0.2, colour = "black") +
  geom_errorbar(
    aes(ymin = AVG_Log_Flux, ymax = AVG_Log_Flux),
    width = 0.5, colour = "black") +
  scale_y_continuous(name = expression('Log10 Flux (kg·m'^-2*'s'^-1*')'),
                     expand = expansion(mult = c(0, 0.05)),
                     limits = c(-11, -4)) +
  scale_x_discrete(name = "Treatment") +
  scale_colour_manual(values = colour_list) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

Mass_T2 <- t2_data %>%
  ggplot(aes(x = Condition, y = Mass_Flux_kgms, colour = Condition)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_errorbar(aes(ymin = AVG_Mass_Flux - SD_Mass_Flux,
                    ymax = AVG_Mass_Flux + SD_Mass_Flux),
                width = 0.2, colour = "black") +
  geom_errorbar(
    aes(ymin = AVG_Mass_Flux, ymax = AVG_Mass_Flux),
    width = 0.5, colour = "black") +
  scale_y_continuous(name = expression('Flux (kg·m'^-2*'s'^-1*')'),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(name = "Treatment") +
  scale_colour_manual(values = colour_list) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display all plots
print(Log_T1); print(Log_T2)
print(Mass_T1); print(Mass_T2)


## Normality testing -------------------------------------------------------
# helper function for normality testing ouptuts
normality_test <- function(data, var_name) {
  data %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      W = if(n() >= 3) shapiro.test(.data[[var_name]])$statistic else NA_real_,
      p_value = if(n() >= 3) shapiro.test(.data[[var_name]])$p.value else NA_real_,
      normal = if(n() >= 3 && p_value > 0.05) "Yes" else "No",
      .groups = "drop"
    )
}

shapiro_log_t1 <- normality_test(t1_data, "Log10_Flux")
shapiro_log_t2 <- normality_test(t2_data, "Log10_Flux")
shapiro_mass_t1 <- normality_test(t1_data, "Mass_Flux_kgms")
shapiro_mass_t2 <- normality_test(t2_data, "Mass_Flux_kgms")

cat("\n=== NORMALITY TESTS ===\n")
cat("T1 Log10:\n"); print(shapiro_log_t1)
cat("T2 Log10:\n"); print(shapiro_log_t2)
cat("T1 Mass:\n"); print(shapiro_mass_t1)
cat("T2 Mass:\n"); print(shapiro_mass_t2)


## Statistical analysis ----------------------------------------------------
# helper function for ANOVA, Dunnett, Tukey tests
stats_by_time <- function(data, timepoint) {
  data$Condition <- relevel(data$Condition, ref = "PDMS")
  
  lm_model <- lm(Log10_Flux ~ Condition, data = data)
  
  dunnett <- glht(lm_model, linfct = mcp(Condition = "Dunnett"))
  tukey <- glht(lm_model, linfct = mcp(Condition = "Tukey"))
  
  cat(sprintf("\n=== %s STATS ===\n", timepoint))
  cat("ANOVA:\n"); print(anova(lm_model))
  cat("Dunnett (PDMS ref):\n"); print(summary(dunnett))
  cat("Tukey:\n"); print(summary(tukey))
  
  cld_tukey <- cld(tukey, alpha = 0.05)
  print(cld_tukey)
  
  return(list(lm = lm_model, tukey = tukey, cld = cld_tukey))
}

t1_stats <- stats_by_time(t1_data, "T1")
t2_stats <- stats_by_time(t2_data, "T2")


# Saving plots and stats --------------------------------------------------
write_csv(Data_summary, "Flux_Summary_T1_T2.csv")
write_csv(t1_stats, "Stats.csv")
ggsave("Log_Flux_T1.pdf", Log_T1, width = 6, height = 8, units = "in")
ggsave("Mass_Flux_T1.pdf", Mass_T1, width = 6, height = 5, units = "in")
ggsave("Log_Flux_T2.pdf", Log_T2, width = 6, height = 5, units = "in")
ggsave("Mass_Flux_T2.pdf", Mass_T2, width = 6, height = 5, units = "in")



