## Data formatting and manipulation of the PDMS bulk permeation experiment
## 
## Bulk diffusion is calculated
## Diffusion rate g/s is calculated
## Values are normalized against the average thickness of the appropriate PDMS control
## Values are converted to the appropriate SI unit for diffusion
## Values are converted into a performance metric compared to the PDMS controls
##
## data is provided in .csv form
##
## 08/03/2025
## Evan Kear

## Set working directory were files are located
setwd("###")

## Load packages
library(tidyverse)


## Reading data ---------------------------------------------------------
# Function to read and rename column
read_and_rename <- function(file_path) {
  read.csv(file_path) %>% rename(Mold_um = Thickness_um)
}

# File directory
data_dir <- "Data"

# File names
files <- list(
  raw_data_ek = "Membrane Bulk Diffusion Water - All Data - EK.csv",
  thickness_data_ek = "Membrane Bulk Diffusion Water - All Thickness Data - EK.csv",
  raw_data_mb = "Membrane Bulk Diffusion Water - All Data - MB.csv",
  thickness_data_mb = "Membrane Bulk Diffusion Water - All Thickness Data - MB.csv"
)

# Load data
datasets <- lapply(files, function(f) read_and_rename(file.path(data_dir, f)))

# Assign variables
raw_data_ek <- datasets$raw_data_ek
thickness_data_ek <- datasets$thickness_data_ek
raw_data_mb <- datasets$raw_data_mb
thickness_data_mb <- datasets$thickness_data_mb


## scaling EK data to same area a MB data ----------------------------------
scale_factor <- 1.131 / 1.474

raw_data_ek <- raw_data_ek %>% 
  mutate(Weight_g = Weight_g * scale_factor)


## Adding thickness averages for each mold to data -----------------------------
# Create values for control thicknesses
Control_Thickness <- data.frame("Control", 0, 0, 0)
names(Control_Thickness) <- c("Name", "Mold_um", "mean", "sd")

# Calculate the average PDMS thickness for normalization
# Function to compute thickness stats and PDMS subset
calculate_thickness_stats <- function(data) {
  stats <- data %>%
    group_by(Name, Mold_um) %>%
    summarise(mean = mean(Actual_Thickness_um, na.rm = TRUE), 
              sd = sd(Actual_Thickness_um, na.rm = TRUE),
              .groups = "drop")
  
  PDMS_stats <- data %>%
    filter(Treatment == "PDMS") %>%
    group_by(Name, Mold_um) %>%
    summarise(mean = mean(Actual_Thickness_um, na.rm = TRUE), 
              sd = sd(Actual_Thickness_um, na.rm = TRUE),
              .groups = "drop")
  
  PDMS_stats <- bind_rows(PDMS_stats, Control_Thickness) # control thickness values are added so it doesnt break with the merge later
  
  list(stats = stats, PDMS = PDMS_stats)
}

# calculate thickness statistics
thickness_ek <- calculate_thickness_stats(thickness_data_ek)
thickness_mb <- calculate_thickness_stats(thickness_data_mb)

thickness_stats_EK <- thickness_ek$stats
PDMS_thickness_EK <- thickness_ek$PDMS
thickness_stats_MB <- thickness_mb$stats
PDMS_thickness_MB <- thickness_mb$PDMS

# Function to merge raw data with PDMS thickness
merge_with_thickness <- function(raw_data, PDMS_thickness) {
  merged_data <- raw_data %>%
    left_join(PDMS_thickness, by = "Mold_um") %>% # Mold_um is now used only for naming conventions
    select(-sd) %>%  # Remove standard deviation if not needed
    rename(Thickness_um = mean, Name = Name.x)
  
  return(merged_data)
}


# Merge data
raw_data_ek <- merge_with_thickness(raw_data_ek, PDMS_thickness_EK)
raw_data_ek <- select(raw_data_ek, -Name.y)
raw_data_mb <- merge_with_thickness(raw_data_mb, PDMS_thickness_MB)
raw_data_mb <- select(raw_data_mb, -Name.y)


# Function to merge and remove NA
merge_and_clean <- function(raw_data, thickness_data) {
  merged_data <- raw_data %>%
    left_join(thickness_data, by = c("Name", "Treatment", "Concentration", "Mold_um", "Replicate")) %>%
    drop_na(Actual_Thickness_um)  # Remove rows where thickness measurement failed
  
  return(merged_data)
}

# Merge and clean data
merged_data_ek <- merge_and_clean(raw_data_ek, thickness_data_ek)
merged_data_mb <- merge_and_clean(raw_data_mb, thickness_data_mb)


## Calculate diffusion -----------------------------------------------------
# Extract initial weights where Time_h == 0
initial_ek <- merged_data_ek %>%
  filter(Time_h == 0) %>%
  rename(StartWeight_g = Weight_g) %>%
  select(-Time_h)

initial_mb <- merged_data_mb %>%
  filter(Time_h == 0) %>%
  rename(StartWeight_g = Weight_g) %>%
  select(-Time_h)

# Add initial weights back into dataset
merged_data_ek <- merged_data_ek %>%
  left_join(initial_ek, by = c(
    "Name", "Treatment", "Concentration", "Thickness_um",
    "Mold_um", "Replicate", "Area_cm2", "Actual_Thickness_um"
  ))

merged_data_mb <- merged_data_mb %>%
  left_join(initial_mb, by = c(
    "Name", "Treatment", "Concentration", "Thickness_um",
    "Mold_um", "Replicate", "Area_cm2", "Actual_Thickness_um"
  )) %>%
 mutate(StartWeight_g = ifelse(is.na(StartWeight_g), 54.22915, StartWeight_g))  # Handle missing values

# Compute diffusion values
calculate_diffusion <- function(data) {
  data %>%
    mutate(
      Diffusion_g = StartWeight_g - Weight_g,
      Normalized_Diffusion_g = ifelse(
        is.na((Actual_Thickness_um / Thickness_um) * Diffusion_g), # here, the diffusion is normalized using the true membrane thickness against the expected thickness (calculated from the PDMS controls) 
        Diffusion_g,
        (Actual_Thickness_um / Thickness_um) * Diffusion_g # NA handling just means where PDMS controls are not valid, it isnt normalized
      )
    ) %>%
    arrange(Time_h, Name, Thickness_um, Replicate) %>%
    group_by(Name, Thickness_um, Replicate) %>%
    mutate(
      Diffusion_Rate_gs = (lag(Weight_g) - Weight_g) / 86400,
      Normalized_Diffusion_Rate_gs = ifelse(
        is.na((Actual_Thickness_um / Thickness_um) * Diffusion_Rate_gs),
        Diffusion_Rate_gs,
        (Actual_Thickness_um / Thickness_um) * Diffusion_Rate_gs
      )
    )
}

diffusion_data_ek <- calculate_diffusion(merged_data_ek)
diffusion_data_mb <- calculate_diffusion(merged_data_mb)

## Conversion to mass flux  ----------------------------------------------------
# Calculate Mass Flux, in units kg/m2/s, currently we are in g/s with no accounting of area
diffusion_data_ek <- diffusion_data_ek %>% 
  mutate(Mass_Flux_kgm2s = ((Normalized_Diffusion_Rate_gs/Area_cm2) * 10))
diffusion_data_mb <- diffusion_data_mb %>% 
  mutate(Mass_Flux_kgm2s = ((Normalized_Diffusion_Rate_gs/Area_cm2) * 10))


## Calculate PDMS averages for each time point/thickness------------------------
# Function to clean data and rename thicknesses
clean_data <- function(data) {
  data <- data %>%   
    mutate(
        Days = Time_h / 24,
        Thickness_um = as.character(Mold_um),
        Thickness_um = case_when(
          Thickness_um == "0"   ~ "Control",
          Thickness_um == "100" ~ "100 µm",
          Thickness_um == "200" ~ "200 µm",
          Thickness_um == "300" ~ "300 µm",
          Thickness_um == "500" ~ "500 µm",
          TRUE ~ Thickness_um
      )
    )
}

data_EK <- clean_data(diffusion_data_ek)
data_MB <- clean_data(diffusion_data_mb)

data_MB <- data_MB %>% # removing replicates that failed, but werent caught by earlier NA removal
  filter(!(Name == "PDMS-0" & Thickness_um == "500 µm" & Replicate == "R7")) %>%
  filter(!(Name == "Control-0" & Thickness_um == "Control" & Replicate == "R5")) %>% 
  filter(!(Name == "Carbopol-15" & Thickness_um == "300 µm" & Replicate == "R7"))


write.csv(data_MB, "Diffusion_Data_MB-All.csv")
write.csv(data_EK, "Diffusion_Data_EK-All.csv")

# function to calculate stats of PDMS, replicates, and treatments
calculate_summary_stats_TIME <- function(data) {
  # Compute PDMS averages at each thickness and time point
  PDMS_TIME_AVG <- data %>%
    na.omit() %>%
    filter(Treatment == "PDMS") %>% 
    group_by(Thickness_um, Time_h, Days) %>%
    summarise(
      PDMS_AVG_FLUX = mean(Mass_Flux_kgm2s, na.rm = TRUE),
      PDMS_AVG_BULK = mean(Normalized_Diffusion_g, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute summary statistics per replicate
  summary_TIME_stats <- data %>%
    na.omit() %>%
    group_by(Name, Thickness_um, Treatment, Replicate, Concentration, Time_h, Days) %>%
    summarise(
      AverageFlux = mean(Mass_Flux_kgm2s, na.rm = TRUE),
      AverageDiffusion = mean(Normalized_Diffusion_g, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute treatment-level summary statistics
  summary_stats2 <- summary_TIME_stats %>%
    group_by(Name, Thickness_um, Treatment, Time_h) %>%
    summarise(
      Treatment_Flux_Avg = mean(AverageFlux, na.rm = TRUE),
      Treatment_Flux_SD = sd(AverageFlux, na.rm = TRUE),
      Treatment_Bulk_Avg = mean(AverageDiffusion, na.rm = TRUE),
      Treatment_Bulk_SD = sd(AverageDiffusion, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge with PDMS averages
  summary_TIME_stats <- left_join(summary_TIME_stats, PDMS_TIME_AVG, by = c("Thickness_um", "Time_h"))
  
  # Merge treatment-level statistics back into the corrected dataset
  summary_TIME_stats <- left_join(summary_TIME_stats, summary_stats2, by = c("Name", "Thickness_um", "Time_h"))
  
  summary_TIME_stats <- summary_TIME_stats %>% 
    select(-Treatment.y, -Days.y) %>% 
    rename(Days = Days.x, Treatment = Treatment.x)
  
  return(summary_TIME_stats)
}


summary_TIME_stats_EK <- calculate_summary_stats_TIME(data_EK) # This dataset can be used for graphing line plots, where time IS considered
summary_TIME_stats_MB <- calculate_summary_stats_TIME(data_MB)

write.csv(summary_TIME_stats_EK, "summary_TIME_stats_EK.csv")
write.csv(summary_TIME_stats_MB, "summary_TIME_stats_MB.csv")

# Function to calculate summary stats, grouped by ALL
calculate_summary_stats <- function(data) {
  # Compute PDMS averages at each thickness and time point
  PDMS_AVG <- data %>%
    na.omit() %>%
    filter(Treatment == "PDMS") %>% 
    group_by(Thickness_um) %>%
    summarise(
      PDMS_AVG_FLUX = mean(Mass_Flux_kgm2s, na.rm = TRUE),
      PDMS_AVG_BULK = mean(Normalized_Diffusion_g, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute summary statistics
  summary_stats <- data %>%
    na.omit() %>%
    group_by(Name, Thickness_um, Treatment, Replicate, Concentration) %>%
    summarise(
      AverageFlux = mean(Mass_Flux_kgm2s, na.rm = TRUE),
      AverageDiffusion = mean(Normalized_Diffusion_g, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compute treatment-level summary statistics
  summary_stats2 <- summary_stats %>%
    group_by(Name, Thickness_um, Treatment) %>%
    summarise(
      Treatment_Flux_Avg = mean(AverageFlux, na.rm = TRUE),
      Treatment_Flux_SD = sd(AverageFlux, na.rm = TRUE),
      Treatment_Bulk_Avg = mean(AverageDiffusion, na.rm = TRUE),
      Treatment_Bulk_SD = sd(AverageDiffusion, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge with PDMS averages
  summary_stats <- left_join(summary_stats, PDMS_AVG, by = c("Thickness_um"))
  
  # Merge treatment-level statistics back into the corrected dataset
  summary_stats <- left_join(summary_stats, summary_stats2, by = c("Name", "Thickness_um"))
  
  summary_stats <- summary_stats %>% 
    select(-Treatment.y) %>% 
    rename(Treatment = Treatment.x)
  
  return(summary_stats)
}

summary_stats_EK <- calculate_summary_stats(data_EK) # This dataset can be used for graphing box/whisker plots, where time IS NOT considered
summary_stats_MB <- calculate_summary_stats(data_MB)

write.csv(summary_stats_EK, "summary_stats_EK.csv")
write.csv(summary_stats_MB, "summary_stats_MB.csv")


## Creating performance metric against PDMS value --------------------------
# Function to calculate performance metric
apply_PDMS_correction_TIME <- function(summary_stats) {
  # Apply PDMS correction
  PDMS_corrected_summary_stats <- summary_stats %>%
    mutate(
      Performance_AverageFlux = AverageFlux / PDMS_AVG_FLUX,
      Performance_AverageDiffusion = AverageDiffusion / PDMS_AVG_BULK
    )
  
  # Calculate treatment-level summary statistics
  summary_stats2 <- PDMS_corrected_summary_stats %>%
    group_by(Name, Thickness_um, Treatment, Time_h, Days) %>%
    summarise(
      Performance_Treatment_Flux_Avg = mean(Performance_AverageFlux, na.rm = TRUE),
      Performance_Treatment_Flux_SD = sd(Performance_AverageFlux, na.rm = TRUE),
      Performance_Treatment_Bulk_Avg = mean(Performance_AverageDiffusion, na.rm = TRUE),
      Performance_Treatment_Bulk_SD = sd(Performance_AverageDiffusion, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge treatment-level statistics back into the corrected dataset
  PDMS_corrected_summary_stats <- left_join(
    PDMS_corrected_summary_stats,
    summary_stats2,
    by = c("Name", "Thickness_um", "Treatment", "Time_h", "Days")
  )
  
  return(PDMS_corrected_summary_stats)
}

Performance_summary_TIME_stats_EK <- apply_PDMS_correction_TIME(summary_TIME_stats_EK)
Performance_summary_TIME_stats_MB <- apply_PDMS_correction_TIME(summary_TIME_stats_MB)

write.csv(Performance_summary_TIME_stats_EK, "Performance_summary_TIME_stats_EK.csv")
write.csv(Performance_summary_TIME_stats_MB, "Performance_summary_TIME_stats_MB.csv")


# Function to calculate performance metric
apply_PDMS_correction <- function(summary_stats) {
  # Apply PDMS correction
  PDMS_corrected_summary_stats <- summary_stats %>%
    mutate(
      Performance_AverageFlux = AverageFlux / PDMS_AVG_FLUX,
      Performance_AverageDiffusion = AverageDiffusion / PDMS_AVG_BULK
    )
  
  # Calculate treatment-level summary statistics
  summary_stats2 <- PDMS_corrected_summary_stats %>%
    group_by(Name, Thickness_um, Treatment) %>%
    summarise(
      Performance_Treatment_Flux_Avg = mean(Performance_AverageFlux, na.rm = TRUE),
      Performance_Treatment_Flux_SD = sd(Performance_AverageFlux, na.rm = TRUE),
      Performance_Treatment_Bulk_Avg = mean(Performance_AverageDiffusion, na.rm = TRUE),
      Performance_Treatment_Bulk_SD = sd(Performance_AverageDiffusion, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge treatment-level statistics back into the corrected dataset
  PDMS_corrected_summary_stats <- left_join(
    PDMS_corrected_summary_stats,
    summary_stats2,
    by = c("Name", "Thickness_um", "Treatment")
  )
  
  return(PDMS_corrected_summary_stats)
}

Performance_summary_stats_EK <- apply_PDMS_correction(summary_stats_EK)
Performance_summary_stats_MB <- apply_PDMS_correction(summary_stats_MB)

write.csv(Performance_summary_stats_EK, "Performance_summary_stats_EK.csv")
write.csv(Performance_summary_stats_MB, "Performance_summary_stats_MB.csv")




