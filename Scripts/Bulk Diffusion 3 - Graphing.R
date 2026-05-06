## Data graphing for PDMS CNC bulk permeation experiment 
##
## data is provided in .csv form,
## created by "Analysis - Bulk Diffusion.R"
##
## 11/03/2024
## Evan Kear

# Load packages
library(tidyverse)
library(gridExtra)
library(patchwork)

# Set working directory
setwd("###")

## Read and prep data ----------------------------------------------------
read_data <- function(prefix) {
  summary <- read.csv(paste0("Performance_summary_stats_", prefix, ".csv"))
  time <- read.csv(paste0("Performance_summary_TIME_stats_", prefix, ".csv"))
  
  list(
    summary = summary,
    time = time
  )
}

data_EK <- read_data("EK")
data_MB <- read_data("MB")

# Filter controls once
data_MB$summary <- data_MB$summary %>% filter(Thickness_um != "Control")
data_MB$time <- data_MB$time %>% filter(Thickness_um != "Control")

# Combine all data
data <- bind_rows(data_EK$summary, data_MB$summary)
data_time <- bind_rows(data_EK$time, data_MB$time)

# Prepare PDMS data
PDMS_data <- bind_rows(
  data_MB$summary %>% 
    filter(Treatment == "PDMS") %>% 
    mutate(Source = "MB"),
  data_EK$summary %>% 
    filter(Treatment == "PDMS") %>% 
    mutate(Source = "EK", Thickness_um = str_replace(Thickness_um, "200 µm", "300 µm"))
)

## Set constants ---------------------------------------------------------------
Treatments <- c("CNC", "Pemulen", "Carbopol", "PVP", "Cellulose", "PDMS")
Thicknesses <- c("100 µm", "200 µm", "300 µm", "500 µm")
ticks <- c(0, 1, 2, 3, 4, 5, 6, 7)

name_levels <- factor(c(
  "PDMS-0", "Control-0",
  "Pemulen-5", "Pemulen-10", "Pemulen-15", 
  "Carbopol-5", "Carbopol-10", "Carbopol-15",
  "PVP-5", "PVP-10", "PVP-15",
  "Cellulose-5", "Cellulose-10", "Cellulose-15",
  "CNC-5", "CNC-10", "CNC-15"
))

data_time$Name <- factor(data_time$Name, levels = name_levels)
data$Name <- factor(data$Name, levels = name_levels)

# Color palettes
color_list_thickness_Pemulen <- c("500 µm" = "#C66BF5", "300 µm" = "#E58FF5", "100 µm" = "#F9ADFF", "Control" = "#000000")
color_list_Pemulen <- c("Control-0" = "#000000", "PDMS-0" = "#A32A3A", "Pemulen-5" = "#F9ADFF", "Pemulen-10" = "#E58FF5", "Pemulen-15" = "#C66BF5")

color_list_thickness_Carbopol <- c("500 µm" = "#0C34AB", "300 µm" = "#6881E6", "100 µm" = "#65B0E6", "Control" = "#000000")
color_list_Carbopol <- c("Control-0" = "#000000", "PDMS-0" = "#A32A3A", "Carbopol-5" = "#65B0E6", "Carbopol-10" = "#6881E6", "Carbopol-15" = "#0C34AB")

color_list_thickness_PVP <- c("500 µm" = "#F58026", "300 µm" = "#F5A03C", "100 µm" = "#F4BE54", "Control" = "#000000")
color_list_PVP <- c("Control-0" = "#000000", "PDMS-0" = "#A32A3A", "PVP-5" = "#F4BE54", "PVP-10" = "#F5A03C", "PVP-15" = "#F58026")

color_list_thickness_Cellulose <- c("500 µm" = "#E0AE51", "300 µm" = "#F5D058", "100 µm" = "#F5E56E", "Control" = "#000000")
color_list_Cellulose <- c("Control-0" = "#000000", "PDMS-0" = "#A32A3A", "Cellulose-5" = "#F5E56E", "Cellulose-10" = "#F5D058", "Cellulose-15" = "#E0AE51")

color_list_thickness_CNC <- c("500 µm" = "#507E4D", "200 µm" = "#A0B76A", "100 µm" = "#B7CF7B", "Control" = "#000000")
color_list_CNC <- c("Control-0" = "#000000", "PDMS-0" = "#A32A3A", "CNC-5" = "#B7CF7B", "CNC-10" = "#A0B76A", "CNC-15" = "#507E4D")

color_list_thickness_PDMS <- c("500 µm" = "#A32A3A", "300 µm" = "#DD394F", "100 µm" = "#FA415A", "Control" = "#000000")
color_list_PDMS <- c("PDMS-0" = "#A32A3A", "Control-0" = "#000000")

colour_list_100 <- c("Pemulen-10" = "#C66BF5", "Carbopol-10" = "#0C34AB", "PVP-10" = "#F58026", "Cellulose-10" = "#E0AE51", "CNC-10" = "#507E4D", "PDMS-0" = "#A32A3A")

## Plotting functions ----------------------------------------------------------
theme_base <- theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.15, 0.8),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

theme_large <- theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# Create time-series line plot with error bars
create_line_plot <- function(.data, y_var, title_suffix, ylim, color_list_name, include_controls = TRUE) {
  # Filter data based on whether controls should be included (uses global treatment/thickness vars)
  plot_data <- if (include_controls) {
    .data %>% filter(Treatment %in% c(treatment, "Control", "PDMS"), Thickness_um %in% c(thickness, "Control"))
  } else {
    .data %>% filter(Treatment == treatment, Thickness_um == thickness)
  }
  
  ggplot(plot_data, aes(Days, {{ y_var }}, colour = Name)) +
    geom_line(linewidth = 0.8) +           # Main trend line
    geom_point(size = 1) +                 # Data points
    geom_errorbar(aes(ymin = {{ y_var }} - SD_var, ymax = {{ y_var }} + SD_var),  # SD error bars
                  width = 0.10, linewidth = 0.5) +
    scale_x_continuous(name = "Days", expand = expansion(add = c(0, 0.05)), breaks = ticks) +
    scale_colour_manual(values = get(color_list_name)) +  # Treatment-specific colors
    labs(title = paste(thickness, "PDMS +", treatment, "Membranes")) +
    coord_cartesian(xlim = c(0, 6), ylim = ylim) +
    guides(colour = guide_legend(position = "inside")) +
    theme_base
}

# Create errorbar + jitter plot for comparing treatment means
create_errorbar_plot <- function(.data, x_var, y_var, color_var, ylim, color_values, title = "", theme = theme_large) {
  ggplot(.data, aes({{ x_var }}, {{ y_var }}, colour = {{ color_var }})) +
    geom_errorbar(aes(ymin = {{ y_var }}, ymax = {{ y_var }}), width = 0.5, colour = "black") +  # Mean position bar
    geom_errorbar(aes(ymin = {{ y_var }} - SD_var, ymax = {{ y_var }} + SD_var), width = 0.2, colour = "black") +  # Mean ± SD bars
    geom_jitter(width = 0.2, size = 3) +      # Individual replicates
    scale_y_continuous(name = expression('Relative Mass Flux'), expand = c(0, 0), limits = ylim) +
    scale_colour_manual(values = color_values) +
    labs(title = title) +
    {{ theme }}
}

## Create line plots -----------------------------------------------------
line_plot_bulk_normalized_list <- list()
line_plot_flux_normalized_list <- list()
line_plot_bulk_performance_list <- list()
line_plot_flux_performance_list <- list()

for (treatment in Treatments) {
  for (thickness in Thicknesses) {
    if (any(data_time$Treatment == treatment & data_time$Thickness_um == thickness)) {
      color_list_name <- paste0("color_list_", treatment)
      
      # Bulk normalized
      line_plot_bulk_normalized_list[[paste0(treatment, "_", thickness, "_bulk")]] <- 
        create_line_plot(data_time, Treatment_Bulk_Avg, "bulk", c(0, 0.5), color_list_name)
      
      # Flux normalized  
      line_plot_flux_normalized_list[[paste0(treatment, "_", thickness, "_flux")]] <- 
        create_line_plot(data_time, Treatment_Flux_Avg, "flux", c(0, 0.000005), color_list_name)
      
      # Bulk performance
      line_plot_bulk_performance_list[[paste0(treatment, "_", thickness, "_bulk_perf")]] <- 
        create_line_plot(data_time, Performance_Treatment_Bulk_Avg, "bulk_perf", c(0, 2), 
                         color_list_name, include_controls = FALSE)
      
      # Flux performance
      line_plot_flux_performance_list[[paste0(treatment, "_", thickness, "_flux_perf")]] <- 
        create_line_plot(data_time, Performance_Treatment_Flux_Avg, "flux_perf", c(0, 2), 
                         color_list_name, include_controls = FALSE)
    }
  }
}

## Generate box plots ------------------------------------------------------
box_plot_bulk_normalized_list <- list()
box_plot_flux_normalized_list <- list()
box_plot_bulk_performance_list <- list()
box_plot_flux_performance_list <- list()

for (treatment in Treatments) {
  for (thickness in Thicknesses) {
    if (any(data$Treatment == treatment & data$Thickness_um == thickness)) {
      color_list_name <- paste0("color_list_", treatment)
      color_values <- get(color_list_name)
      
      filt_data <- data %>% filter(Treatment %in% c(treatment, "PDMS"), Thickness_um == thickness)
      
      # Bulk normalized
      box_plot_bulk_normalized_list[[paste0(treatment, "_", thickness, "_bulk")]] <- 
        create_errorbar_plot(filt_data, Name, Treatment_Bulk_Avg, Name, c(0, 0.3), color_values)
      
      # Flux normalized
      box_plot_flux_normalized_list[[paste0(treatment, "_", thickness, "_flux")]] <- 
        create_errorbar_plot(filt_data, Name, Treatment_Flux_Avg, Name, c(0, 0.000005), color_values)
      
      # Bulk performance
      box_plot_bulk_performance_list[[paste0(treatment, "_", thickness, "_bulk_perf")]] <- 
        create_errorbar_plot(filt_data, Name, Performance_Treatment_Bulk_Avg, Name, c(0, 2), 
                             color_values, theme = theme_base)
      
      # Flux performance
      box_plot_flux_performance_list[[paste0(treatment, "_", thickness, "_flux_perf")]] <- 
        create_errorbar_plot(filt_data, Name, Performance_Treatment_Flux_Avg, Name, c(0, 2), 
                             color_values, theme = theme_large)
    }
  }
}

## Thickness plots ---------------------------------------------------------
box_plot_thickness_flux_performance_list <- list()

for (treatment in Treatments) {
  color_values <- get(paste0("color_list_thickness_", treatment))
  
  for (conc in unique(data$Concentration[data$Treatment == treatment])) {
    if (any(data$Treatment == treatment & data$Concentration == conc)) {
      filt_data <- data %>% 
        filter(Treatment %in% c(treatment, "PDMS"), Concentration == conc)
      
      box_plot_thickness_flux_performance_list[[paste0(treatment, "_", conc, "_flux")]] <- 
        create_errorbar_plot(filt_data, Thickness_um, Performance_AverageFlux, Thickness_um, 
                             c(0.5, 2.5), color_values)
    }
  }
}

# Special PDMS normalized box plot
PDMS_normalized_box <- data_MB$summary %>%
  filter(Treatment == "PDMS") %>%
  create_errorbar_plot(Thickness_um, AverageFlux, Thickness_um, c(0, 0.000005), 
                       color_list_thickness_PDMS, theme = theme_large) +
  scale_y_continuous(name = expression('Mass Flux (kgs' ^ '-1' * 'm' ^ '-2' * ')'))

# All excipients plots (preserving exact filtering logic)
create_all_excipients_plot <- function(thickness_filter) {
  filtered_data <- data %>%
    filter((Thickness_um %in% thickness_filter) & (Concentration %in% c(10, 0))) %>%
    filter(Treatment %in% c("PDMS", "CNC", "Cellulose", "PVP", "Pemulen", "Carbopol")) %>%  # Removes controls explicitly
    mutate(Treatment = factor(Treatment, levels = c("PDMS", "CNC", "Cellulose", "PVP", "Pemulen", "Carbopol")))
  
  colour_list <- c("PDMS" = "#A32A3A", "CNC" = "#507E4D", "Cellulose" = "#E0AE51", 
                   "PVP" = "#F58026", "Pemulen" = "#C66BF5", "Carbopol" = "#0C34AB")
  
  create_errorbar_plot(filtered_data, Treatment, Performance_AverageFlux, Treatment, 
                       c(0, 2.2), colour_list, 
                       title = paste0("Thickness: ", paste(thickness_filter, collapse = "/")), 
                       theme_large + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}

plot_100 <- create_all_excipients_plot("100 µm")
plot_300 <- create_all_excipients_plot(c("200 µm", "300 µm"))
plot_500 <- create_all_excipients_plot("500 µm")

## Save example plots (preserving specific saves) --------------------------
output_dir <- ("###")
ggsave("Cellulose_Conc_flux.pdf", box_plot_flux_performance_list$`Cellulose_300 µm_flux_perf`, width = 6, height = 5, units = "in")
ggsave("Carbopol_10_flux.pdf", box_plot_thickness_flux_performance_list$Carbopol_10_flux, width = 6, height = 5, units = "in")
ggsave("PDMS_10_flux.pdf", box_plot_thickness_flux_performance_list$PDMS_0_flux, width = 6, height = 5, units = "in")  # Using generated PDMS
ggsave("PDMS_Thickness_Normalized_Flux.pdf", PDMS_normalized_box, width = 6, height = 5)
ggsave("Relative_Flux_100.pdf", plot_100, height = 5, width = 6)
