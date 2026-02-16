# Gathering data from supporting material published with
# Knock et al., (2021), ‘Key Epidemiological Drivers and Impact of Interventions in the 2020 SARS-CoV-2 Epidemic in England’.
# Science Translational Medicine, 16(602), eabg4262

# Age disaggregated data relating to severity of SARS-Cov-2 in 2020
# age groups converted to single year ages ranging from 0 to 100

# Data extracted
# (1) IFR and IHR by age (see Table s9 and Figure 4C)
# (2) Proportion asymptomatic/symptomatic by age
# (3) Care pathway probabilities


################################################################################
# (1) IFR and IHR (see Table s9 and Figure 4C)
################################################################################

# Data copied from Table S9
ihr_ifr_by_age <- data.frame( age_group = c( rep("0_to_5",5), rep("5_to_10",5),rep("10_to_15",5),rep("15_to_20",5),rep("20_to_25",5),rep("25_to_30",5),rep("30_to_35",5),rep("35_to_40",5),rep("40_to_45",5),rep("45_to_50",5),rep("50_to_55",5),rep("55_to_60",5),rep("60_to_65",5),rep("65_to_70",5),rep("70_to_75",5),rep("75_to_80",5),rep("80+",21),"Care_home_residents")
                              ,age_single_year = c(seq(0,100,1),"Care_home_residents")
                              , IHR = c(rep(0.23,5),rep(0.25,5),rep(0.077,5),rep(0.033,5),rep(0.057,5),rep(0.29,5),rep(0.40,5),rep(0.57,5),rep(0.92,5),rep(1.4,5),rep(1.8,5),rep(2.9,5),rep(6.7,5),rep(7.9,5),rep(17,5),rep(30,5),rep(30,21),19)
                              , IHR_lower = c(rep(0.18,5),rep(0.20,5),rep(0.061,5),rep(0.026,5),rep(0.045,5),rep(0.23,5),rep(0.31,5),rep(0.45,5),rep(0.72,5),rep(1.1,5),rep(1.4,5),rep(2.3,5),rep(5.3,5),rep(6.2,5),rep(13,5),rep(23,5),rep(24,21),4.8)
                              , IHR_upper =c(rep(0.29,5),rep(0.32,5),rep(0.097,5),rep(0.041,5),rep(0.072,5),rep(0.37,5),rep(0.50,5),rep(0.72,5),rep(1.2,5),rep(1.8,5),rep(2.3,5),rep(3.7,5),rep(8.5,5),rep(9.9,5),rep(21,5),rep(37,5),rep(38,21),35)
                              , IFR = c(rep(0.009,5),rep(0.010,5),rep(0.003,5),rep(0.001,5),rep(0.002,5),rep(0.013,5),rep(0.020,5),rep(0.031,5),rep(0.057,5),rep(0.106,5),rep(0.165,5),rep(0.331,5),rep(0.939,5),rep(1.344,5),rep(3.359,5),rep(6.761,5),rep(7.622,21),23.133)
                              , IFR_lower =c(rep(0.006,5),rep(0.007,5),rep(0.002,5),rep(0.001,5),rep(0.002,5),rep(0.009,5),rep(0.014,5),rep(0.022,5),rep(0.040,5),rep(0.076,5),rep(0.118,5),rep(0.242,5),rep(0.693,5),rep(0.996,5),rep(2.524,5),rep(5.047,5),rep(5.722,21),14.434)
                              , IFR_upper =c(rep(0.014,5),rep(0.015,5),rep(0.005,5),rep(0.002,5),rep(0.004,5),rep(0.019,5),rep(0.028,5),rep(0.043,5),rep(0.078,5),rep(0.145,5),rep(0.222,5),rep(0.441,5),rep(1.247,5),rep(1.759,5),rep(4.396,5),rep(8.841,5),rep(10.024,21),14.434)
                              )

library(ggplot2)
library(tidyr)
library(dplyr)

### First option

# Convert age_single_year to a factor to handle "Care_home_residents"
ihr_ifr_by_age$age_single_year <- factor( ihr_ifr_by_age$age_single_year, levels = c(as.character(seq(0,100,1)), "Care_home_residents"))

# Reshape to long format
ihr_ifr_long <- pivot_longer( ihr_ifr_by_age, cols = c("IHR", "IFR"), names_to = "Measure", values_to = "Value")

# Plot
ggplot( ihr_ifr_long, aes(x = age_single_year, y = Value, color = Measure, group = Measure)) +
  geom_line() +
  labs(x = "Age", y = "Rate (%)", title = "IHR and IFR by Age") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### Second option
# Remove "Care_home_residents"
ihr_ifr_by_age_filtered <- ihr_ifr_by_age %>%
  filter( age_single_year != "Care_home_residents" )

# Ensure age_single_year is a factor with numeric levels
ihr_ifr_by_age_filtered$age_single_year <- as.numeric( as.character( ihr_ifr_by_age_filtered$age_single_year ) )

# Reshape to long format
ihr_ifr_long <- pivot_longer( ihr_ifr_by_age_filtered, cols = c("IHR", "IFR"), names_to = "Measure", values_to = "Value")

# Create x-axis breaks and labels for every 5 years
breaks_5yr <- seq(0, 100, 5)

# Plot
ggplot( ihr_ifr_long, aes(x = age_single_year, y = Value, color = Measure, group = Measure)) +
  geom_line(linewidth = 1.5) +
  scale_x_continuous(
    breaks = breaks_5yr,
    labels = as.character(breaks_5yr)
  ) +
  labs(x = "Age (years)", y = "Rate(%)", title = "") +
  theme_minimal() +
  theme(  axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12)
        , axis.title.x = element_text(size = 12)
        , axis.text.y = element_text( size = 12)
        , axis.title.y = element_text(size = 12)
        , plot.title = element_text(size = 16)
        , legend.title = element_text(size = 16)
        , legend.text = element_text(size = 16)
        )

### Option with CIs
ihr_ifr_by_age_filtered <- ihr_ifr_by_age %>%
  filter( age_single_year != "Care_home_residents" ) %>%
  mutate(
    IHR_lower = IHR_lower#IHR * 0.9,
    ,IHR_upper = IHR_upper#IHR * 1.1,
    ,IFR_lower = IFR_lower#IFR * 0.9,
    ,IFR_upper = IFR_upper#IFR * 1.1
  )

# Convert age_single_year to numeric
ihr_ifr_by_age_filtered$age_single_year <- as.numeric( as.character( ihr_ifr_by_age_filtered$age_single_year ) )

# Reshape to long format for plotting
ihr_ifr_long <- ihr_ifr_by_age_filtered %>%
  pivot_longer(
    cols = c(IHR, IFR),
    names_to = "Measure",
    values_to = "Value"
  ) %>%
  mutate(
    lower = ifelse(Measure == "IHR", IHR_lower, IFR_lower),
    upper = ifelse(Measure == "IHR", IHR_upper, IFR_upper)
  )

# Define breaks every 5 years
breaks_5yr <- seq(0, 100, 5)

# Plot with confidence intervals, increased font size and line width
ggplot( ihr_ifr_long, aes(x = age_single_year, y = Value, color = Measure, fill = Measure, group = Measure)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_x_continuous(
    breaks = breaks_5yr,
    labels = as.character(breaks_5yr)
  ) +
  labs(x = "Age (years)", y = "Rate (%)", title = "Age-stratified estimates of disease severity") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)
  )

saveRDS( ihr_ifr_by_age
         , "~/mscape/1_epidemic_modelling/Knock et al 2021/KD_outputs_for_nbpmscape/knock_ihr_ifr_by_age.rds")

################################################################################
# (2) Proportion asymptomatic/symptomatic by age
################################################################################
# See section 2.3.4 in Supplementary Material of Knock et al (2021)
# Note: not including value at upper end of range
symptomatic_prob_by_age <- data.frame( age_group = c( rep("0_to_5",5), rep("5_to_10",5),rep("10_to_15",5),rep("15_to_20",5),rep("20_to_25",5),rep("25_to_30",5),rep("30_to_35",5),rep("35_to_40",5),rep("40_to_45",5),rep("45_to_50",5),rep("50_to_55",5),rep("55_to_60",5),rep("60_to_65",5),rep("65_to_70",5),rep("70_to_75",5),rep("75_to_80",5),rep("80+",21),"Care_home_residents", "Care_home_workers")
                                     , age_single_year = c(seq(0,100,1),"Care_home_residents","Care_home_workers")
                                     , symptomatic_prob = c( rep(0.25,5), rep(0.26875,5),rep(0.325,5),rep(0.41875,5),rep(0.55,83))
                                     )

saveRDS( symptomatic_prob_by_age
         , "~/mscape/1_epidemic_modelling/Knock et al 2021/KD_outputs_for_nbpmscape/knock_symptomatic_prob_by_age.rds")

################################################################################
# (3) Care pathway probabilities
################################################################################

# Load data from Knock et al (2021) Supplementary file found at https://www.science.org/doi/10.1126/scitranslmed.abg4262#supplementary-materials

care_pathway_prob_by_age_group <- read.csv("~/mscape/1_epidemic_modelling/Knock et al 2021/abg4262_data_file_s1.csv")
care_pathway_prob_by_age_group <- as.data.frame( t( care_pathway_prob_by_age_group ) )
colnames( care_pathway_prob_by_age_group ) <- care_pathway_prob_by_age_group[1,]
care_pathway_prob_by_age_group <- care_pathway_prob_by_age_group[-1,]

# Add column for age group
care_pathway_prob_by_age_group <- cbind( care_pathway_prob_by_age_group
                                   , age_group = rep( NA, nrow(care_pathway_prob_by_age_group) ) )

for(i in 1: nrow( care_pathway_prob_by_age_group)){
  ag_old_1 <- rownames( care_pathway_prob_by_age_group )[i]
  ag_old_2 <- gsub("^X","", ag_old_1 )
  ag_new   <- gsub("\\.","_", ag_old_2 )
  care_pathway_prob_by_age_group$age_group[i] <- ag_new
}

# Convert to single year age df
care_pathway_prob_by_age <- as.data.frame( matrix( data = NA
                                                   , nrow = 101
                                                   , ncol = ncol( care_pathway_prob_by_age_group ) + 1
                                                   ) )
colnames( care_pathway_prob_by_age )[1:ncol(care_pathway_prob_by_age_group)] = colnames(care_pathway_prob_by_age_group)
colnames( care_pathway_prob_by_age )[ncol(care_pathway_prob_by_age)] = "age_single_year"

# Fill single year age rows with data from age group rows
for ( j in 1:nrow( care_pathway_prob_by_age_group )){

  # Identify age range in age group
  ag_lower <- as.numeric( strsplit( care_pathway_prob_by_age_group$age_group[j], "_")[[1]][1] )
  ag_upper <- as.numeric( strsplit( care_pathway_prob_by_age_group$age_group[j], "_")[[1]][3] )

  # Fill in data
  care_pathway_prob_by_age[ ag_lower : ag_upper+1, 1:ncol(care_pathway_prob_by_age_group)] <- care_pathway_prob_by_age_group[j,]

  # Add single year age to row
  care_pathway_prob_by_age$age_single_year[ ag_lower : ag_upper+1 ] <- seq( ag_lower, ag_upper, 1 )

}

# Rearrange columns
care_pathway_prob_by_age <- care_pathway_prob_by_age[,c(11,10,seq(1,9))]

saveRDS( care_pathway_prob_by_age
         , "~/mscape/1_epidemic_modelling/Knock et al 2021/KD_outputs_for_nbpmscape/knock_care_pathway_prob_by_age.rds")
