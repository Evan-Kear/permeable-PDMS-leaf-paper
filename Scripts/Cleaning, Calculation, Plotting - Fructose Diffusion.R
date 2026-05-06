## This is a script to clean and graph platereader data from 
## fructose diffusion experiments, with readings for OD and YFP.
## Data is inputted in RAW form, cleaned into long format and exported.
## Curves are visually inspected for all samples/treatments to ensure data robustness. 
## The standard curve information is extracted, averaged and plotted to 
## create the needed curves and equations, based on AUC. 
## The same is done for the respective experiment data. 
##
## Data is provided in .csv form
## 17/11/2025
## Evan Kear

## Define working directory
wd <- ("XXX")

## Packages
library(tidyverse)
library(DescTools)
library(patchwork)

## Set current experiment for analysis
excipient <- "PVP"

## Set the dilution ratio
dilution <- 1/4

## Read membrane thicknesses
setwd(wd)
thickness_data <- read.csv("Fructose Diffusion Membrane Thickness.csv")

## Set new working directory for the excipient
setwd(paste(wd, excipient, "/", sep = ""))

## Read RAW data
YFP_RAW <- read.csv(paste("YFP_Raw.csv", sep =""))

## Cleaning Data -----------------------------------------------------------
## Function to clean data. (AI generated bassed on self written code - Perplexity)
process_plate_reader_data <- function(df) {
  require(tidyverse)
  
  # Pivot longer and convert types
  df_long <- df %>%
    pivot_longer(
      cols = 2:ncol(df),
      names_to = "Time",
      values_to = "YFP"
    ) %>%
    mutate(
      Time = as.integer(str_extract(Time, "\\\\d+")),  # More robust numeric extraction
      YFP = as.numeric(YFP)
    )
  
  # Seperate main data
  all_data <- df_long %>% 
    separate(X, c("Condition", "Sample Point", "Replicate", "Technical_Replicate"), sep = "_")
  
  main_data <- all_data %>%   
    filter(Condition == excipient | Condition == 'Control' | Condition == "Donor")
  
  # Extract standards
  standards_data <- all_data %>% 
    filter(Condition == "STD") %>% 
    select(-Replicate) %>% 
    rename(Concentration = `Sample Point`)
  
  # Calculate the 0 mM values at each timepoint, these need to be removed from all observations
  STD_0mM <- standards_data %>% 
    filter(Concentration == "0.0 mM") %>%
    group_by(Time) %>% 
    summarise(STD_0 = mean(YFP))
  
  # Subtract blanks and adjust for 0 mM signal in main data
  # (at this point I realise that the media blank is accounted for in the 
  #  0 mM standard, so I dont need to subtract both)
  adj_data <- main_data %>% 
    left_join(STD_0mM, by = c("Time")) %>% 
    mutate(YFP_ADJ = (YFP - STD_0)) %>% 
    select(-YFP, -STD_0)
  
  donor_data <- adj_data %>% 
    filter(Condition == 'Donor')
  
  # Subtract blanks and adjust for 0 mM signal in STD data
  # (at this point I realise that the media blank is accounted for in the 
  #  0 mM standard, so I dont need to subtract both)
  adj_STD <- standards_data %>% 
    left_join(STD_0mM, by = c("Time")) %>% 
    mutate(YFP_ADJ = (YFP - STD_0)) %>% 
    select(-YFP, -STD_0) %>% 
    filter(Condition == "STD")
  
  # Write cleaned data to CSVs
  data_file_name <- paste(excipient, "YFP Data Clean.csv")
  std_file_name <- paste(excipient, "YFP Standards Clean.csv")
  
  write.csv(adj_data, data_file_name, row.names = FALSE)
  write.csv(adj_STD, std_file_name, row.names = FALSE)
  
  # Return dataset as a list
  return(list(adj_data = adj_data, adj_STD = adj_STD))
}

# Run function (df)
YFP <- process_plate_reader_data(YFP_RAW)
Standards <- YFP$adj_STD
Main_data <- YFP$adj_data

# Replicate lists
replicates <- unique(Main_data$Replicate) 
sample_replicates <- replicates
sample_replicates <- replicates[-((length(replicates) - 2):length(replicates))]
control_replicates <- tail(unique(Main_data$Replicate), 3)

## Visually Checking the STDs ----------------------------------------------
concentrations <- unique(Standards$Concentration)

color_list <- c("6.25 mM" = "#FF2445", "3.125 mM" = "#FF3256", "1.56 mM" = "#FF4167", 
                "0.78 mM" = "#FF4F78", "0.39 mM" = "#FF5E89", "0.2 mM" = "#FF6C9B", 
                "0.1 mM" = "#FF7AAC", "0.05 mM" = "#FF89BD", "0.025 mM" = "#FF97CE", 
                "0.0125 mM" = "#FFA6DF","0.00625 mM" = "#FFB4F0", "0.0 mM" = "#4f474a")

# Set levels
Standards$Concentration <- factor(
  Standards$Concentration,
  levels =
    c("6.25 mM",
      "3.125 mM",
      "1.56 mM",
      "0.78 mM",
      "0.39 mM",
      "0.2 mM",
      "0.1 mM",
      "0.05 mM",
      "0.025 mM",
      "0.0125 mM",
      "0.00625 mM",
      "0.0 mM")
)

plot_std_replicate <- function(std_data, tech_rep, plot_title, y_limit) {
  std_data %>% 
    filter(Technical_Replicate == tech_rep) %>% 
    ggplot(aes(Time, YFP_ADJ, colour = as.factor(Concentration), group = Concentration)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = y_limit, xlim = c(0, 480)) + 
    labs(
      title = plot_title,
      colour = "Concentration"
    ) +
    theme_classic() +
    scale_colour_manual(values = color_list) +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
}

STD_TR1_Plot <- plot_std_replicate(Standards, "Tr1", "TR1 Concentration Curves", c(0, 8000))
STD_TR1_Plot

STD_TR2_Plot <- plot_std_replicate(Standards, "Tr2", "TR2 Concentration Curves", c(0, 8000))
STD_TR2_Plot

## Merging STD values ------------------------------------------------------
# If the STD pass visual checks, can average them. 

Standards <- Standards %>% 
  group_by(Time, Concentration) %>% 
  mutate(YFP_ADJ = mean(YFP_ADJ)) %>% 
  filter(Technical_Replicate == "Tr1") %>%
  select(-Technical_Replicate)

STD_Plot <- Standards %>% 
  ggplot(aes(Time, YFP_ADJ, colour = as.factor(Concentration), group = Concentration)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 6000), xlim = c(0, 480)) + 
  labs(
    title = "Standard Concentration Curves",
    colour = "Concentration"
  ) +
  theme_classic() +
  scale_colour_manual(values = color_list) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
STD_Plot

ggsave(paste0(excipient,"_Standard_Curves.pdf"), STD_Plot, width = 8, height = 6)

## Plotting T1 technical replicates of each chamber ------------------------
# Making lists and then drawing plots for the tech reps
T1_replicates_tests <- list()
T1_replicates_controls <- list()

plot_technical_replicates <- function(data, replicate_values, condition_value, sample_point_value, output_list_name) {
  plot_list <- list()
  
  for (rep in replicate_values) {
    p <- data %>% 
      filter(Replicate == rep) %>% 
      filter(Condition == condition_value) %>% 
      filter(`Sample Point` == sample_point_value) %>% 
      ggplot(aes(Time, YFP_ADJ, colour = Technical_Replicate)) +
      geom_line() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = c(0, 6000), xlim = c(0, 480)) + 
      labs(
        title = paste0("Tech Replicates for: ", rep),
      ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)
      )
    
    plot_list[[paste0("Replicate_", rep)]] <- p
  }
  
  assign(output_list_name, plot_list, envir = parent.frame())
}

plot_technical_replicates(Main_data, sample_replicates, excipient, "T1", "T1_replicates_tests")
T1_replicates_tests$Replicate_R9

wrap_plots(T1_replicates_tests)

# Use the following to remove bad replicates. 
T1_Test_Data <- Main_data %>% 
  filter(Condition == paste0(excipient), `Sample Point` == "T1", Time <= 480,)

# Calculate median and deviation
T1_Test_Data <- T1_Test_Data %>%
  group_by(Replicate, Time) %>%
  mutate(
    median_YFP = median(YFP_ADJ, na.rm = TRUE),
    deviation = abs(YFP_ADJ - median_YFP)
  )

# Thresholding and flag abnormal replicates
threshold <- 3500 # Adjust threshold based on data
T1_Test_Data <- T1_Test_Data %>%
  mutate(is_normal = deviation < threshold)

# Filter to keep only normal
T1_Test_Data_Clean <- T1_Test_Data %>% filter(is_normal)

# Averaging the normal replicates
T1_Test_Data_Final <- T1_Test_Data_Clean %>% 
  group_by(Time, Replicate) %>% 
  mutate(YFP = mean(YFP_ADJ)) %>% 
  filter(Technical_Replicate == "Tr1") %>% 
  select(-YFP_ADJ, -Technical_Replicate, -median_YFP, -deviation, -is_normal)

## Plotting T2 technical replicates of each chamber ------------------------
# Making lists and then drawing plots for the tech reps
T2_replicates_tests <- list()
T2_replicates_controls <- list()
replicates <- unique(Main_data$Replicate)

plot_technical_replicates(Main_data, sample_replicates, excipient, "T2", "T2_replicates_tests")
T2_replicates_tests$Replicate_R7

wrap_plots(T2_replicates_tests)

# Use the following to remove bad replicates. 
T2_Test_Data <- Main_data %>% 
  filter(Condition == paste0(excipient), `Sample Point` == "T2", Time <= 480,
         Replicate != "R7")

# Calculate median and deviation
T2_Test_Data <- T2_Test_Data %>%
  group_by(Replicate, Time) %>%
  mutate(
    median_YFP = median(YFP_ADJ, na.rm = TRUE),
    deviation = abs(YFP_ADJ - median_YFP)
  )

# Thresholding and flag abnormal replicates
threshold <- 3500 # Adjust threshold based on data
T2_Test_Data <- T2_Test_Data %>%
  mutate(is_normal = deviation < threshold)

# Filter to keep only normal
T2_Test_Data_Clean <- T2_Test_Data %>% filter(is_normal)

# Averaging the normal replicates
T2_Test_Data_Final <- T2_Test_Data_Clean %>% 
  group_by(Time, Replicate) %>% 
  mutate(YFP = mean(YFP_ADJ)) %>% 
  filter(Technical_Replicate == "Tr1") %>% 
  select(-YFP_ADJ, -Technical_Replicate, -median_YFP, -deviation, -is_normal)

## Plotting T1 control replicates of each chamber ------------------------
# Making lists and then drawing plots for the tech reps
T1_replicates_controls <- list()

plot_technical_replicates(Main_data, control_replicates, "Control", "T1", "T1_replicates_controls")
T1_replicates_controls$Replicate_R3

wrap_plots(T1_replicates_controls)

# Use the following to remove bad replicates. 
T1_Control_Data <- Main_data %>% 
  filter(Condition == "Control", `Sample Point` == "T1", Time <= 480,
         Replicate == "R9")

# Calculate median and deviation
T1_Control_Data <- T1_Control_Data %>%
  group_by(Replicate, Time) %>%
  mutate(
    median_YFP = median(YFP_ADJ, na.rm = TRUE),
    deviation = abs(YFP_ADJ - median_YFP)
  )

# Thresholding and flag abnormal replicates
threshold <- 3500 # Adjust threshold based on data
T1_Control_Data <- T1_Control_Data %>%
  mutate(is_normal = deviation < threshold)

# Filter to keep only normal
T1_Control_Data_Clean <- T1_Control_Data %>% filter(is_normal)

# Averaging the normal replicates
T1_Control_Data_Final <- T1_Control_Data_Clean %>% 
  group_by(Time, Replicate) %>% 
  mutate(YFP = mean(YFP_ADJ)) %>% 
  filter(Technical_Replicate == "Tr1") %>% 
  select(-YFP_ADJ, -Technical_Replicate, -median_YFP, -deviation, -is_normal)

# Plotting T2 control replicates of each chamber ------------------------
# Making lists and then drawing plots for the tech reps
T2_replicates_controls <- list()
replicates <- unique(Main_data$Replicate)

plot_technical_replicates(Main_data, control_replicates, "Control", "T2", "T2_replicates_controls")
T2_replicates_controls$Replicate_R3

wrap_plots(T2_replicates_controls)

# Use the following to remove bad replicates. 
T2_Control_Data <- Main_data %>% 
  filter(Condition == "Control", `Sample Point` == "T2", Time <= 480, 
         Replicate == "R9")

# Calculate median and deviation
T2_Control_Data <- T2_Control_Data %>%
  group_by(Replicate, Time) %>%
  mutate(
    median_YFP = median(YFP_ADJ, na.rm = TRUE),
    deviation = abs(YFP_ADJ - median_YFP)
  )

# Thresholding and flag abnormal replicates
threshold <- 3500 # Adjust threshold based on data
T2_Control_Data <- T2_Control_Data %>%
  mutate(is_normal = deviation < threshold)

# Filter to keep only normal
T2_Control_Data_Clean <- T2_Control_Data %>% filter(is_normal)

# Averaging the normal replicates
T2_Control_Data_Final <- T2_Control_Data_Clean %>% 
  group_by(Time, Replicate) %>% 
  mutate(YFP = mean(YFP_ADJ)) %>% 
  filter(Technical_Replicate == "Tr1") %>% 
  select(-YFP_ADJ, -Technical_Replicate, -median_YFP, -deviation, -is_normal)

## Merging Cleaned Data and plotting ---------------------------------------
# Merging cleaned tech rep data 

Clean_Data <- bind_rows(
  T1_Test_Data_Final,
  # T2_Test_Data_Final,
  #T1_Control_Data_Final,
  # T2_Control_Data_Final
)

plot_curve_lines <- function(data, sample_point_value, condition_value, plot_title) {
  data %>% 
    filter(`Sample Point` == sample_point_value) %>% 
    filter(Condition == condition_value) %>% 
    ggplot(aes(Time, YFP, colour = as.factor(Replicate), group = Replicate)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 7000), xlim = c(0, 480)) + 
    labs(
      title = plot_title,
      colour = "Concentration"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
}

Curves_Plot_T1 <- plot_curve_lines(Clean_Data, "T1", excipient, "T1 Curves")
Curves_Plot_T1

Curves_Plot_Control_T1 <- plot_curve_lines(Clean_Data, "T1", "Control", "T1 Curves")
Curves_Plot_Control_T1

Curves_Plot_T2 <- plot_curve_lines(Clean_Data, "T2", excipient, "T2 Curves")
Curves_Plot_T2

Curves_Plot_Control_T2 <- plot_curve_lines(Clean_Data, "T2", "Control", "T2 Curves")
Curves_Plot_Control_T2

Curve_plots <- list(Curves_Plot_T1, Curves_Plot_T2, Curves_Plot_Control_T1, Curves_Plot_Control_T2)
wrap_plots(Curve_plots)

ggsave(paste0(excipient,"_Sample_Plots.pdf"), Curve_plots, width = 8, height = 6)

## AUC Calculations For Standards --------------------------------------------------------
# calculating AUC AVERGAE 
Standards <- Standards %>%
  mutate(Treatment = recode(Concentration,
                            "6.25 mM" = "6.25",
                            "3.125 mM" = "3.125",
                            "1.56 mM" = "1.56",
                            "0.78 mM" = "0.78",
                            "0.39 mM" = "0.39",
                            "0.2 mM" = "0.2",
                            "0.1 mM" = "0.1",
                            "0.05 mM" = "0.05",
                            "0.025 mM" = "0.025",
                            "0.0125 mM" = "0.0125",
                            "0.00625 mM" = "0.00625",
                            "0.0 mM" = "0.0"
  ))

AUC_Standards_All <- Standards %>%
  filter(Time <= 480 & Treatment != 0) %>% 
  group_by(Treatment) %>%
  summarize(
    AUC_YFP = AUC(x = Time, y = pmax(YFP_ADJ, 0), method = "trapezoid")
  )

AUC_Standards_All$Treatment <- as.numeric(as.character(AUC_Standards_All$Treatment))

write.csv(AUC_Standards_All, paste0(excipient, "_Standard_AUC_Values.csv"))
# Here, determine what range is usable for the curve

## Fitting monod and curve -------------------------------------------------
# Fit data to monod equation
monod_fit_all <- AUC_Standards_All %>%
  { nls(AUC_YFP ~ (Vmax * Treatment) / (Ks + Treatment),
        data = .,
        start = list(Vmax = max(.$AUC_YFP), Ks = median(.$Treatment))) }

# Create new data for predictions (fine grid for smooth curve)
newdata <- tibble(Treatment = seq(min(AUC_Standards_All$Treatment), max(AUC_Standards_All$Treatment), length.out = 100))

# Predict AUC_YFP using the fitted Monod model
newdata <- newdata %>%
  mutate(AUC_YFP = predict(monod_fit_all, newdata))

# Plot points and Monod fitted curve
plot_auc_Standards_All <- AUC_Standards_All %>% 
  ggplot(aes(x = Treatment, y = AUC_YFP)) +
  geom_point(size = 3) +
  geom_line(data = newdata, aes(x = Treatment, y = AUC_YFP), color = "blue", size = 1) +
  scale_y_continuous(name = "AUC", expand = c(0, 0)) +
  scale_x_continuous(name = "96 Well Concentration [mM]", expand = c(0, 0)) +
  labs(title = paste0("MRE612 Standard Curve ", excipient)) +
  coord_cartesian(xlim = c(0, 1.56), ylim = c(0, 2000000)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

plot_auc_Standards_All

ggsave(paste0(excipient, "_Standards_AUC_MONOD.png"), plot_auc_Standards_All, width = 8, height = 7, units = "in")

## AUC Calculations For Standards LOWER RANGE --------------------------------------------------------
# calculating AUC AVERGAE 

AUC_Standards_LOWRANGE <- Standards %>%
  filter(Time <= 480 & Treatment != 0.0125 & Treatment != 0.025 & Treatment != 0.00625) %>% 
  group_by(Treatment) %>%
  summarize(
    AUC_YFP = AUC(x = Time, y = pmax(YFP_ADJ, 0), method = "trapezoid")
  )

AUC_Standards_LOWRANGE$Treatment <- as.numeric(as.character(AUC_Standards_LOWRANGE$Treatment))

write.csv(AUC_Standards_LOWRANGE, paste0(excipient, "_Standard_AUC_Values_LOWRANGE.csv"))


## Fitting monod and curve LOWER RANGE -------------------------------------------------
# Fit data to monod equation
monod_fit_LOWRANGE <- AUC_Standards_LOWRANGE %>%
  { nls(AUC_YFP ~ (Vmax * Treatment) / (Ks + Treatment),
        data = .,
        start = list(Vmax = max(.$AUC_YFP), Ks = median(.$Treatment))) }

# Create new data for predictions (fine grid for smooth curve)
newdata <- tibble(Treatment = seq(min(AUC_Standards_LOWRANGE$Treatment), max(AUC_Standards_LOWRANGE$Treatment), length.out = 100))

# Predict AUC_YFP using the fitted Monod model
newdata <- newdata %>%
  mutate(AUC_YFP = predict(monod_fit_LOWRANGE, newdata))

# Plot points and Monod fitted curve
plot_auc_Standards_LOWRANGE <- AUC_Standards_LOWRANGE %>% 
  ggplot(aes(x = Treatment, y = AUC_YFP)) +
  geom_point(size = 3) +
  geom_line(data = newdata, aes(x = Treatment, y = AUC_YFP), color = "blue", size = 1) +
  scale_y_continuous(name = "AUC", expand = c(0, 0)) +
  scale_x_continuous(name = "96 Well Concentration [mM]", expand = c(0, 0)) +
  labs(title = paste0("MRE612 Standard Curve ", excipient)) +
  coord_cartesian(xlim = c(0, 1.56), ylim = c(0, 2000000)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

plot_auc_Standards_LOWRANGE

ggsave(paste0(excipient, "_Standards_AUC_MONOD_LOWRANGE.png"), plot_auc_Standards_LOWRANGE, width = 8, height = 7, units = "in")


## Checking equations against predicted ------------------------------------
# Coefficients (Vmax and Ks) from the monod model
params_all <- coef(monod_fit_all)
Vmax_all <- params_all["Vmax"]
Ks_all   <- params_all["Ks"]

# Coefficients (Vmax and Ks) from the monod model
params_LOW <- coef(monod_fit_LOWRANGE)
Vmax_LOW <- params_LOW["Vmax"]
Ks_LOW   <- params_LOW["Ks"]

# Create data frame and save to CSV
params_df <- data.frame(
  Model = c("All", "LowRange"),
  Vmax = c(Vmax_all, Vmax_LOW),
  Ks = c(Ks_all, Ks_LOW)
)

write_csv(params_df, paste0(excipient, "_monod_parameters.csv"))

# Check standards and save
AUC_Standards_LOWRANGE_CHECKS <- AUC_Standards_LOWRANGE %>% 
  mutate(Prediction = (AUC_YFP * Ks_LOW) / (Vmax_LOW - AUC_YFP)) %>% 
  mutate(Chamber_Conc = Prediction/4) %>% 
  mutate(True_Chamber_Conc = Treatment/4)

AUC_Standards_CHECKS <- AUC_Standards_All %>% 
  mutate(Prediction = (AUC_YFP * Ks_all) / (Vmax_all - AUC_YFP)) %>% 
  mutate(Chamber_Conc = Prediction/4) %>% 
  mutate(True_Chamber_Conc = Treatment/4)

# This allows the appropriate model to be selected for the prediction

write.csv(AUC_Standards_LOWRANGE_CHECKS, paste0(excipient, "AUC_Standards_LOWRANGE_CHECKS.csv"))
write.csv(AUC_Standards_CHECKS, paste0(excipient, "AUC_Standards_CHECKS.csv"))

## Calculating Sample AUC --------------------------------------------------
Sample_AUC <- Clean_Data %>%
  filter(Time <= 480) %>% 
  group_by(Condition, Replicate, `Sample Point`) %>%
  summarize(
    AUC_YFP = AUC(x = Time, y = pmax(YFP, 0), method = "trapezoid")
  )

# Add prediction to the AUC data
Sample_AUC <- Sample_AUC %>%
  mutate(Prediction = (AUC_YFP * Ks_LOW) / (Vmax_LOW - AUC_YFP))

# Here is a good time to check the AUC outputs, 
# do they make sense?
# The equation may force values into negatives, which isnt right.
# Are the controls similar values, does this indicate leaking?
# Are the values similar?


## Converting total diffusion to rate --------------------------------------
# Firstly, the data collected is still in its diluted form, 
# so the prediction needs to be adjusted back to the chamber concentration
Sample_AUC <- Sample_AUC %>% 
  mutate(Chamber_Conc = Prediction/dilution)

write.csv(Sample_AUC, paste0(excipient, "_Sample_Concentration_Prediction.csv"))

# this converts concentration of 1ml of fructose solution into grams
Diffusion_metrics <- Sample_AUC %>% 
  mutate(Grams = Chamber_Conc * 0.00018016) %>% 
  filter(Chamber_Conc >= 0)

# Due to the variation in membrane thickness, an actual and target thickness
# need to be established.
target_thickness = mean(thickness_data$Thickness)

replicate_thickness = thickness_data %>% 
  filter(Membrane == excipient) %>% 
  select(-Membrane)

# I set the control thickness to be equal to the target so no changes will be made
Diffusion_metrics <- Diffusion_metrics %>% 
  left_join(replicate_thickness) %>%
  mutate(Thickness = if_else(Condition == "Control", target_thickness, Thickness))

Diffusion_metrics <- Diffusion_metrics %>% 
  mutate(
    Normalized_Diffusion_Rate_gs = (Thickness / target_thickness * Grams) / 28800,  # Diffusion Rate calculation
    Mass_Flux_kgms = (Normalized_Diffusion_Rate_gs / 1.130973355) * 10,  # Area/unit conversion
    Log10_Flux = log10(Mass_Flux_kgms)
  )

write.csv(Diffusion_metrics, paste0(excipient, "_Diffusion_metrics.csv"))

Diffusion_metrics_summary <- Diffusion_metrics %>% 
  drop_na() %>% 
  group_by(Condition, `Sample Point`) %>% 
  summarise(
    AVG_Log_FLux = mean(Log10_Flux),
    SD_Log_FLux = sd(Log10_Flux)
  )

write.csv(Diffusion_metrics_summary, paste0(excipient, "_Diffusion_metrics_summary.csv"))