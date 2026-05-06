## This script is for analyzing the CFU counts from the membrane growth experiments
##
## Evan Kear
## 28/10/2025

setwd("XXX")

# Libraries
library(tidyverse)

# Read data
all_data <- read.csv("All Survival Data.csv")
all_data <- all_data %>% mutate(Humidity = ifelse(Humidity == "H20", 100, Humidity))

# Split into experiments
unique_experiments <- unique(all_data$Experiment)

# Empty list to hold data frames
experiment_data_list <- list()

# Unique experiments and filter
for (i in seq_along(unique_experiments)) {
  experiment_num <- unique_experiments[i]
  experiment_data_list[[i]] <- all_data %>% filter(Experiment == experiment_num)
  # Optionally, set names for easy access
  names(experiment_data_list)[i] <- paste0("data_", experiment_num)
}


# Long Term Experiment ----------------------------------------------------
# Clean data for long term experiment
data_longterm <- experiment_data_list[["data_5"]] %>% 
  mutate(Count = ifelse(is.na(as.numeric(Count)), NA, as.numeric(Count))) %>% 
  drop_na()

data_longterm <- data_longterm %>% 
  mutate(Log_CFU_ml_cm = log10(as.numeric(CFU...ml...cm) + 1))


# Calculating avg and SD

data_long_average <- data_longterm %>% 
  group_by(Time, Humidity, Condition) %>% 
  summarise(Avg = mean(Log_CFU_ml_cm),
            SD = sd(Log_CFU_ml_cm))

# Plotting

longterm_survival_plot <- data_long_average %>% 
  ggplot(aes(Time, Avg, colour = Condition, fill = Condition)) +
  geom_line() +
  geom_ribbon(aes(ymin = Avg - SD, ymax = Avg + SD), alpha = 0.2, colour = NA) +  # add envelope
  scale_y_continuous(name = expression(log[10](CFU/ml/cm)), limits = c(0,5.5), expand = c(0, 0)) +
  scale_x_continuous(name = "Time (h)", expand = c(0, 0)) +
  labs(title = "Membrane bacterial survival 100% RH") +
  
  theme_minimal()

longterm_survival_plot



# Medium Term Survival ----------------------------------------------------
# Clean data for long term experiment
data_medium <- experiment_data_list[["data_4"]] %>% 
  mutate(Count = ifelse(is.na(as.numeric(Count)), NA, as.numeric(Count))) %>% 
  drop_na()

data_medium <- data_medium %>% 
  mutate(Log_CFU_ml_cm = log10(as.numeric(CFU...ml...cm) + 1))


# Calculating avg and SD

data_medium_average <- data_medium %>% 
  group_by(Time, Humidity, Condition) %>% 
  summarise(Avg = mean(Log_CFU_ml_cm),
            SD = sd(Log_CFU_ml_cm))

# Plotting

mediumterm_survival_plot <- data_medium_average %>% 
  ggplot(aes(Time, Avg, colour = Condition, fill = Condition)) +
  geom_line() +
  geom_ribbon(aes(ymin = Avg - SD, ymax = Avg + SD), alpha = 0.2, colour = NA) +  # add envelope
  scale_y_continuous(name = expression(log[10](CFU/ml/cm)), limits = c(0,5.5), expand = c(0, 0)) +
  scale_x_continuous(name = "Time (h)", expand = c(0, 0)) +
  labs(title = "Membrane bacterial survival 100% RH") +
  
  theme_minimal()

mediumterm_survival_plot

# Humidity Survival ----------------------------------------------------
# Clean data for long term experiment
data_Humidity <- experiment_data_list[["data_2"]] %>% 
  mutate(Count = ifelse(is.na(as.numeric(Count)), NA, as.numeric(Count))) %>% 
  drop_na()

data_Humidity <- data_Humidity %>% 
  mutate(Log_CFU_ml_cm = log10(as.numeric(CFU...ml...cm) + 1))


# Calculating avg and SD

data_Humidity_average <- data_Humidity %>% 
  group_by(Time, Humidity, Condition) %>% 
  summarise(Avg = mean(Log_CFU_ml_cm),
            SD = sd(Log_CFU_ml_cm))

# Plotting

Humidity_survival_plot <- data_Humidity_average %>% 
  filter(Condition == "Growth") %>% 
  ggplot(aes(Time, Avg, colour = as.factor(Humidity))) +
  geom_line() +
  geom_ribbon(aes(ymin = Avg - SD, ymax = Avg + SD), alpha = 0.2, colour = NA) +  # add envelope
  scale_y_continuous(name = expression(log[10](CFU/ml/cm)), limits = c(0,5.5), expand = c(0, 0)) +
  scale_x_continuous(name = "Time (h)", expand = c(0, 0)) +
  labs(title = "Membrane bacterial survival 100% RH") +
  
  theme_minimal()

Humidity_survival_plot
