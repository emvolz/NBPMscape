# Description: 
#   Script to assess the effect of increasing symptomatic proportion on 
#   time to detection using sampling in different healthcare settings
# Date: October 2026

# Using
# RStudio 2023.06.1+524 "Mountain Hydrangea" Release (547dcf861cac0253a8abb52c135e44e02ba407a1, 2023-07-07) for windows
# Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) RStudio/2023.06.1+524 Chrome/110.0.5481.208 Electron/23.3.0 Safari/537.36

R.version.string
#[1] "R version 4.5.0 (2025-04-11 ucrt)"

# ── Packages ──────────────────────────────────────────────────────────
install.packages(c("car", "dunn.test", "rstatix", "ggplot2", "tidyr"))
install.packages("rstudioapi")
library(car); library(dunn.test); library(rstatix); library(ggplot2); library(tidyr)
library(rstudioapi)

packageVersion("car") # 3.1.5
packageVersion("dunn.test") # 1.4.2
packageVersion("rstatix") # 1.1.0
packageVersion("ggplot2") # 4.0.3
packageVersion("tidyr") # 1.3.2
packageVersion("rstudioapi") # 0.19.0

# ── Prepare Data: Read and manipulate ──────────────────────────

# Read data computed in 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/6_varying_symptomatic_fraction_.jl'
base_folder <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

td_symp_df <- read.csv( paste0( base_folder, "/covid_like_varied_sympt_prop/td_symp_df_GP_ICU.csv" ) )

# Filter data by scenario 1 (300 samples per week in winter)
colnames(td_symp_df)
#[1] "symptomatic_proportion" "sampling_method"        "scenario_number"        "scenario_name"         
#[5] "sim_rep"                "TD" 
unique(td_symp_df$scenario_number) #[1] 13  1 25 73 75 74
# Scenario numbers represent seasons (winter/summer), different numbers of samples per week
# and also different numbers of ICU sites
#                                     300 spw   500 spw     1000 spw
#_____________________________________________________________________
# Primary care (RCGP) / ICU 8 sites |   1        13          25
# ICU 26 sites                      |  73        74          75
#
# See 'NBPMscape/config/config file record.xlsx' for additional scenario descriptions


# Rename sampling_method to show the number of ICU sites
for (i in 1:nrow(td_symp_df)){ #i=1
  snum    = td_symp_df[i,"scenario_number"]
  smethod = td_symp_df[i,"sampling_method"]
  if ( snum %in% c(1,13,25) & smethod == "ICU"){
    td_symp_df[i,"sampling_method"] = "ICU_sites_8"
  } else if ( snum %in% c(73,74,75) & smethod == "ICU"){
    td_symp_df[i,"sampling_method"] = "ICU_sites_26"
  }
}

# Check successful
unique(td_symp_df$sampling_method) # [1] "RCGP" "ICU_sites_26" "ICU_sites_8" 
unique(td_symp_df$scenario_number) # [1] 13  1 25 73 75 74

View( subset(td_symp_df, scenario_number %in% c(73,74,75) & sampling_method == "ICU" ) )
# Should produce empty df

# Set up dataframe to hold differences and ratios
td_symp_diff_ratio_df <- data.frame(
   symptomatic_proportion = character() #numeric()
  ,sampling_method_1      = character()
  ,sampling_method_2      = character()
  ,scenario_number        = numeric()
  ,scenario_name          = character()
  ,sim_rep                = numeric()
  ,td_diff                = numeric()
  ,td_ratio               = numeric()
  ,log_td_ratio           = numeric()
)

symp_props = unique(td_symp_df$symptomatic_proportion) #[1] "base" "base_minus_20pp" "base_plus_20pp"  "base_plus_40pp"  "all_100" 
unique(td_symp_df$sampling_method) # [1] "RCGP" "ICU_sites_26" "ICU_sites_8"
samp_methods_list = list(c("RCGP","ICU_sites_8"),c("RCGP","ICU_sites_26"))
sim_reps_total = 1000
scenario_pairs_list = list(c(1,73),c(13,74),c(25,75))
scenario_numbers = c(1,13,25)
scenario_numbers_icu26 = c(73,74,75)

# Function to compute difference, ratio and log(ratio), or return NA
safe_diff_ratio <- function(a, b) {
  # Check if both are numeric and not NA/NaN
  if (is.numeric(a) && is.numeric(b) && !is.na(a) && !is.na(b) && !is.nan(a) && !is.nan(b)) {
    return( c(a-b, a/b, log(a/b)) )
  } else {
    return(NA_real_)  # Return numeric NA
  }
}

counter = 0
for (sp in symp_props){ #sp = symp_props[2]
  for (j in 1:length(scenario_numbers)){ # j=1 
    for (sm in samp_methods_list){ # sm = samp_methods_list[[1]] # sm[1] # sm[2]
      for (i in 1:sim_reps_total){ # i=1
        
        counter = counter +1
        #print(counter)
        
        scen_number_1 <- ifelse( sm[1] == "ICU_sites_26", scenario_numbers_icu26[j], scenario_numbers[j] )
        scen_number_2 <- ifelse( sm[2] == "ICU_sites_26", scenario_numbers_icu26[j], scenario_numbers[j] ) #scenario_numbers[j] # Assumes that ICU_sites_26 scenario number is only ever the numerator in ratio and first number in the diff
        
        # Subset TD results
        results_sub_1 <- subset(  td_symp_df
                                  , sampling_method == sm[1] & 
                                    symptomatic_proportion == sp &
                                    scenario_number == scen_number_1 &
                                    sim_rep == i )
        results_sub_2 <- subset(  td_symp_df
                                  , sampling_method == sm[2] & 
                                    symptomatic_proportion == sp & 
                                    scenario_number == scen_number_2 & 
                                    sim_rep == i )
        # Compute difference and ratio
        safe_diff_ratio_results <- safe_diff_ratio( results_sub_1$TD, results_sub_2$TD )
        td_diff      <- safe_diff_ratio_results[1]
        td_ratio     <- safe_diff_ratio_results[2]
        log_td_ratio <- safe_diff_ratio_results[3] 
        # Create new row for diff and ratio results df
        temp_row = data.frame(   symptomatic_proportion = sp
                               , sampling_method_1 = sm[1]
                               , sampling_method_2 = sm[2]
                               , scenario_number = scen_number_1
                               , sim_rep = i
                               , td_diff = td_diff
                               , td_ratio = td_ratio
                               , log_td_ratio = log_td_ratio
                               )
        # Add row to df
        td_symp_diff_ratio_df <- rbind(td_symp_diff_ratio_df, temp_row)
      }
    }
  }
}

# Check df compilation is correct
test_1 = subset(  td_symp_df
                  , sampling_method == "RCGP" & 
                    symptomatic_proportion == "base" &
                    scenario_number == 1 &
                    sim_rep == 1 )
print(test_1)
#       symptomatic_proportion  sampling_method   scenario_number scenario_name             sim_rep TD
# 1001  base                    RCGP              1               gp_tds_winter_300samples  1       42.36945

test_2 = subset(  td_symp_df
                  , sampling_method == "ICU_sites_26" & 
                    symptomatic_proportion == "base" &
                    scenario_number == 73 &
                    sim_rep == 1 )
print(test_2)
#       symptomatic_proportion  sampling_method scenario_number scenario_name                                                sim_rep   TD
# 3001  base                    ICU_sites_26    73              icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites  1         45.65841

test_1$TD / test_2$TD == subset( td_symp_diff_ratio_df, symptomatic_proportion == "base" &
                                   sampling_method_1 == "RCGP" &
                                   sampling_method_2 == "ICU_sites_26" &
                                   scenario_number == 1 &
                                   sim_rep == 1 )$td_ratio
# TRUE
print(test_1$TD / test_2$TD)
#[1] 0.9279659

test_1$TD - test_2$TD == subset( td_symp_diff_ratio_df, symptomatic_proportion == "base" &
                                   sampling_method_1 == "RCGP" &
                                   sampling_method_2 == "ICU_sites_26" &
                                   scenario_number == 1 &
                                   sim_rep == 1 )$td_diff
# TRUE
print(test_1$TD - test_2$TD)
#[1] -3.288964

# Add text to describe the sampling methods used in diff and ratio values
td_symp_diff_ratio_df$sampling_methods = paste0(td_symp_diff_ratio_df$sampling_method_1
                                                ,"_v_"
                                                ,td_symp_diff_ratio_df$sampling_method_2)

# Add alternative description of the symptomatic proportion for plotting labels
for ( i in 1:nrow(td_symp_diff_ratio_df)){
  if ( td_symp_diff_ratio_df$symptomatic_proportion[i] == "base_minus_20pp"){
    td_symp_diff_ratio_df$symp_prop[i] <- -20 }
  else if (td_symp_diff_ratio_df$symptomatic_proportion[i] == "base"){
    td_symp_diff_ratio_df$symp_prop[i] <- 0 }
  else if (td_symp_diff_ratio_df$symptomatic_proportion[i] == "base_plus_20pp"){
    td_symp_diff_ratio_df$symp_prop[i] <- 20}
  else if ( td_symp_diff_ratio_df$symptomatic_proportion[i] == "base_plus_40pp"){
    td_symp_diff_ratio_df$symp_prop[i] <- 40}
  else if ( td_symp_diff_ratio_df$symptomatic_proportion[i] == "all_100"){
    td_symp_diff_ratio_df$symp_prop[i] <- 100}
  
}


## Compute means and medians
for (sm in c("RCGP","ICU_sites_26")){
  for (sn in c(1,13,25)){
    for (sp in c("base_minus_20pp","base","base_plus_20pp","base_plus_40pp","all_100")){
      # ICU sites 26 uses different scenario numbers for the same scenario
      if (sm == "ICU_sites_26"){
        if (sn == 1 ){sn = 73}
        if (sn == 13 ){sn = 74}
        if (sn == 25 ){sn = 75}
        }
      temp_td_df <- subset( td_symp_df, sampling_method == sm &
                                        scenario_number == sn &
                                        symptomatic_proportion == sp &
                                        !is.infinite(TD) )
      print( paste0( sm," ",sn," ",sp," mean = ", round(mean(temp_td_df$TD),2), ", median = ", round(median(temp_td_df$TD),2) ) )
    }
  }
}

#[1] "RCGP          1 base_minus_20pp   mean = 55.93, median = 56.2"
#[1] "RCGP          1 base              mean = 49.98, median = 50.48"
#[1] "RCGP          1 base_plus_20pp    mean = 44.56, median = 44.64"
#[1] "RCGP          1 base_plus_40pp    mean = 40.56, median = 40.7"
#[1] "RCGP          1 all_100           mean = 38.29, median = 38.53"
#[1] "RCGP          13 base_minus_20pp  mean = 53.51, median = 53.44"
#[1] "RCGP          13 base             mean = 47.22, median = 47.3"
#[1] "RCGP          13 base_plus_20pp   mean = 42.29, median = 42.47"
#[1] "RCGP          13 base_plus_40pp   mean = 38.36, median = 38.5"
#[1] "RCGP          13 all_100          mean = 36,    median = 36.3"
#[1] "RCGP          25 base_minus_20pp  mean = 49.22, median = 49.35"
#[1] "RCGP          25 base             mean = 43.48, median = 44.09"
#[1] "RCGP          25 base_plus_20pp   mean = 38.71, median = 38.86"
#[1] "RCGP          25 base_plus_40pp   mean = 35.39, median = 35.33"
#[1] "RCGP          25 all_100          mean = 33.79, median = 33.83"
#[1] "ICU_sites_26  73 base_minus_20pp  mean = 50.42, median = 50.63"
#[1] "ICU_sites_26  73 base             mean = 44.73, median = 44.66"
#[1] "ICU_sites_26  73 base_plus_20pp   mean = 41.01, median = 41.09"
#[1] "ICU_sites_26  73 base_plus_40pp   mean = 37.31, median = 37.32"
#[1] "ICU_sites_26  73 all_100          mean = 36.14, median = 36.23"
#[1] "ICU_sites_26  74 base_minus_20pp  mean = 46.79, median = 47.17"
#[1] "ICU_sites_26  74 base             mean = 42.11, median = 42.39"
#[1] "ICU_sites_26  74 base_plus_20pp   mean = 38.09, median = 38.12"
#[1] "ICU_sites_26  74 base_plus_40pp   mean = 35.04, median = 35.28"
#[1] "ICU_sites_26  74 all_100          mean = 33.79, median = 34.1"
#[1] "ICU_sites_26  75 base_minus_20pp  mean = 44.93, median = 45.03"
#[1] "ICU_sites_26  75 base             mean = 39.82, median = 40.17"
#[1] "ICU_sites_26  75 base_plus_20pp   mean = 36.17, median = 36.19"
#[1] "ICU_sites_26  75 base_plus_40pp   mean = 33.43, median = 33.52"
#[1] "ICU_sites_26  75 all_100          mean = 32.27, median = 32.53"

# Compute mean and median differences
results_list <- list()  # container to collect each row of results
row_counter <- 1

for (sm in c("RCGP_v_ICU_sites_26")){ # sm = "RCGP_v_ICU_sites_26"
  for (sn in c(1,13,25)){ # sn = 1
    for (sp in c("base_minus_20pp","base","base_plus_20pp","base_plus_40pp","all_100")){ #sp = "base_minus_20pp"
      temp_diff_ratio_df <- subset( td_symp_diff_ratio_df, sampling_methods == sm &
                                      scenario_number == sn &
                                      symptomatic_proportion == sp &
                                      !is.infinite(td_diff) &
                                      !is.nan(td_diff) )
      
      # compile mean and median stats into a single row instead of printing
      results_list[[row_counter]] <- data.frame(
        sampling_method          = sm,
        scenario_number          = sn,
        symptomatic_proportion   = sp,
        mean_diff                = mean( temp_diff_ratio_df$td_diff ),
        mean_ratio               = mean( temp_diff_ratio_df$td_ratio ),
        mean_log_td_ratio        = mean( temp_diff_ratio_df$log_td_ratio ),
        median_diff              = median( temp_diff_ratio_df$td_diff ),
        median_ratio             = median( temp_diff_ratio_df$td_ratio ),
        median_log_td_ratio      = median( temp_diff_ratio_df$log_td_ratio ),
        stringsAsFactors = FALSE
      )
      
      row_counter <- row_counter + 1
    }
  }
}

# combine all rows into one dataframe
td_diff_ratio_summary_df <- do.call(rbind, results_list)

# view the compiled results
print(td_diff_ratio_summary_df)

# save to CSV
write.csv(td_diff_ratio_summary_df
          , paste0( base_folder,"/covid_like_varied_sympt_prop/td_diff_ratio_summary.csv")
          , row.names = FALSE)


### Focus on RCGP vs ICU_sites_26 sampling
## Split out sampling scenarios
R_v_I26_scenario_1 <- subset( td_symp_diff_ratio_df, sampling_methods == "RCGP_v_ICU_sites_26" &
                                                   scenario_number == 1 &
                                                   !is.infinite(td_diff)
                            )
unique(R_v_I26_scenario_1$symptomatic_proportion)
unique(R_v_I26_scenario_1$sampling_method_1)
unique(R_v_I26_scenario_1$sampling_method_2)
unique(R_v_I26_scenario_1$scenario_number)
unique(R_v_I26_scenario_1$sampling_methods)
unique(R_v_I26_scenario_1$symp_prop)


R_v_I26_scenario_13 <- subset( td_symp_diff_ratio_df, sampling_methods == "RCGP_v_ICU_sites_26" &
                              scenario_number == 13 &
                               !is.infinite(td_diff)
)

R_v_I26_scenario_25 <- subset( td_symp_diff_ratio_df, sampling_methods == "RCGP_v_ICU_sites_26" &
                              scenario_number == 25 &
                               !is.infinite(td_diff)
)

# ════════════════════════════════════════════════════════════════════════
# STEP 1: CHECK NORMALITY (Shapiro-Wilk per condition)
# ════════════════════════════════════════════════════════════════════════
# If p > 0.05 → normality assumption holds → proceed with ANOVA
# If not then normality assumption does not hold and so use the Kruskal-Wallis test.
# This tests whether the rank distributions differ across groups (ie symptomatic proportion).
library(dplyr)

# RCGP vs ICU (26 sites) for scenario 1 (winter, 300 samples per week)
cat("\n── Normality: Difference ──\n")
R_v_I26_scenario_1 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(td_diff)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion     p_value
#<chr>                        <dbl>
#1 all_100                0.0000169  
#2 base                   0.00000574 
#3 base_minus_20pp        0.0232     
#4 base_plus_20pp         0.000000247
#5 base_plus_40pp         0.0000770

# Conclusion: Requirement for normality not met for any symptomatic proportions

cat("\n── Normality: Ratio ──\n")
R_v_I26_scenario_1 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(td_ratio)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion  p_value
#<chr>                     <dbl>
#1 all_100                2.41e-24
#2 base                   4.85e-28
#3 base_minus_20pp        6.46e-33
#4 base_plus_20pp         2.41e-41
#5 base_plus_40pp         1.52e-37

# Conclusion: Requirement for normality not met for any symptomatic proportions

#** Log(ratio) more likely to be normally distributed **
cat("\n── Normality: Log(Ratio) ──\n")
R_v_I26_scenario_1 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(log_td_ratio)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion  p_value
#<chr>                     <dbl>
#1 all_100                2.08e-13
#2 base                   7.59e-16
#3 base_minus_20pp        6.25e-16
#4 base_plus_20pp         8.62e-18
#5 base_plus_40pp         1.46e-15

# Conclusion: Requirement for normality not met for any symptomatic proportions


# RCGP vs ICU (26 sites) for scenario 13 (winter, 500 samples per week)
cat("\n── Normality: Difference ──\n")
R_v_I26_scenario_13 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(td_diff)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion    p_value
#<chr>                       <dbl>
#1 all_100                0.0000891
#2 base                   0.000365 
#3 base_minus_20pp        0.0000255
#4 base_plus_20pp         0.0000222
#5 base_plus_40pp         0.0000347

# Conclusion: Requirement for normality (at 5% significance level) not 
#             met for any symptomatic proportions


cat("\n── Normality: Ratio ──\n")
R_v_I26_scenario_13 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(td_ratio)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion  p_value
#<chr>                     <dbl>
#1 all_100                2.74e-28
#2 base                   1.41e-27
#3 base_minus_20pp        7.89e-37
#4 base_plus_20pp         1.98e-29
#5 base_plus_40pp         1.12e-29

# Conclusion: Requirement for normality not met for any symptomatic proportions

cat("\n── Normality: Log(Ratio) ──\n")
R_v_I26_scenario_13 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(log_td_ratio)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion  p_value
#<chr>                     <dbl>
#1 all_100                2.83e-15
#2 base                   3.25e-13
#3 base_minus_20pp        1.37e-18
#4 base_plus_20pp         1.28e-14
#5 base_plus_40pp         2.05e-14

# Conclusion: Requirement for normality not met for any symptomatic proportions

# RCGP vs ICU (26 sites) for scenario 25 (winter, 1000 samples per week)

cat("\n── Normality: Difference ──\n")
R_v_I26_scenario_25 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(td_diff)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion     p_value
#<chr>                        <dbl>
#1 all_100                0.000938  
#2 base                   0.0000158 
#3 base_minus_20pp        0.00000276
#4 base_plus_20pp         0.0000195 
#5 base_plus_40pp         0.00000303

# Conclusion: Requirement for normality (at 5% significance level) not
#             met for any symptomatic proportions


cat("\n── Normality: Log(Ratio) ──\n")
R_v_I26_scenario_25 |> group_by(symptomatic_proportion) |>
  summarise(p_value = shapiro.test(log_td_ratio)$p.value) |> print()
# A tibble: 5 × 2
#symptomatic_proportion  p_value
#<chr>                     <dbl>
#1 all_100                1.13e-13
#2 base                   8.58e-20
#3 base_minus_20pp        1.77e-18
#4 base_plus_20pp         1.79e-19
#5 base_plus_40pp         9.12e-16

# Conclusion: Requirement for normality not met for any symptomatic proportions


# ════════════════════════════════════════════════════════════════════════
# STEP 2: CHECK HOMOGENEITY OF VARIANCE (Levene's Test)
# ════════════════════════════════════════════════════════════════════════
# p > 0.05 → equal variances → standard ANOVA is valid
# If p-value is below level of significance then variances are NOT equal and 
# ANOVA should not be used

cat("\n── Levene's Test: Difference ──\n")
leveneTest(td_diff ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_1)
#── Levene's Test: Difference ──
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
#group    4  31.616 < 2.2e-16 ***
#      4970                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conclusion: p-value below significance level and so there are significant 
#             differences between the variances of differences between TDs at 
#             different symptomatic proportion

cat("\n── Levene's Test: Log(Ratio) ──\n")
leveneTest(log_td_ratio ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_1)
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value Pr(>F)
#group    4  0.7028 0.5899
#      4970 

# Conclusion: p-value IS ABOVE 5% significance level and so there are NOT significant 
#             differences between the variances of the log_TD_ratios between TDs at 
#             different symptomatic proportions


cat("\n── Levene's Test: Difference ──\n")
leveneTest(td_diff ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_13)
#Levene's Test for Homogeneity of Variance (center = median)
#         Df F value    Pr(>F)    
#group    4  19.618 4.908e-16 ***
#       4984
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conclusion: p-value below significance level and so there are significant 
#             differences between the variances of differences between TDs at 
#             different symptomatic proportion

cat("\n── Levene's Test: Log(Ratio) ──\n")
leveneTest(log_td_ratio ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_13)
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value Pr(>F)
#group    4  1.3716 0.2411
#      4984                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conclusion: p-value IS ABOVE 5% significance level and so there are NOT significant 
#             differences between the variances of the log_TD_ratios between TDs at 
#             different symptomatic proportions


cat("\n── Levene's Test: Difference ──\n")
leveneTest(td_diff ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_25)
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
#group    4  19.751 3.798e-16 ***
#     4989                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conclusion: p-value below significance level and so there are significant 
#             differences between the variances of differences between TDs at 
#             different symptomatic proportion

cat("\n── Levene's Test: Log(Ratio) ──\n")
leveneTest(log_td_ratio ~ as.factor(symptomatic_proportion), data = R_v_I26_scenario_25)
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value  Pr(>F)  
#group    4  2.4761 0.0422 *
#       4989                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conclusion: p-value below significance level and so there are significant 
#             differences between the variances of differences between TDs at 
#             different symptomatic proportion

#### ANOVA is not appropriate  ####
# because not all of the data is normally distributed and 
# not all of the variances are the same

# ════════════════════════════════════════════════════════════════════════
# STEP 3A: ONE-WAY ANOVA — DIFFERENCE
# ════════════════════════════════════════════════════════════════════════
#cat("\n── ANOVA: Difference ──\n")
#anova_diff <- aov(difference ~ condition, data = data)
#summary(anova_diff)

# ════════════════════════════════════════════════════════════════════════
# STEP 3B: ONE-WAY ANOVA — LOG(RATIO)
# ════════════════════════════════════════════════════════════════════════
#cat("\n── ANOVA: Log(Ratio) ──\n")
#anova_ratio <- aov(log_ratio ~ condition, data = data)
#summary(anova_ratio)

# ════════════════════════════════════════════════════════════════════════
# STEP 4: POST-HOC — TUKEY HSD (if ANOVA is significant)
# ════════════════════════════════════════════════════════════════════════
#cat("\n── Tukey HSD: Difference ──\n")
#TukeyHSD(anova_diff)

#cat("\n── Tukey HSD: Log(Ratio) ──\n")
#TukeyHSD(anova_ratio)

# ════════════════════════════════════════════════════════════════════════
# STEP 5: NON-PARAMETRIC ALTERNATIVE — KRUSKAL-WALLIS
# ════════════════════════════════════════════════════════════════════════
# If Kruskal-Wallis test is below level of significance then cannot reject
# the null hypothesis that the medians of all groups are equal


# RCGP vs ICU (26 sites) for scenario 1 (winter, 300 samples per week)

cat("\n── Kruskal-Wallis: Difference ──\n")
kruskal.test(td_diff ~ symptomatic_proportion, data = R_v_I26_scenario_1)
#Kruskal-Wallis rank sum test
#data:  td_diff by symptomatic_proportion
#Kruskal-Wallis chi-squared = 91.978, df = 4, p-value < 2.2e-16

cat("\n── Kruskal-Wallis: Ratio (raw) ──\n")
kruskal.test(td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_1)
#Kruskal-Wallis rank sum test
#data:  td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 34.009, df = 4, p-value = 7.42e-07

cat("\n── Kruskal-Wallis: log(Ratio) ──\n")
kruskal.test(log_td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_1)
#Kruskal-Wallis rank sum test
#data:  log_td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 34.009, df = 4, p-value = 7.42e-07

# Conclusion: Kruskal-Wallis p-value < 0.05 for td_diff, td_ratio and log_td_ratio
#             therefore CANNOT reject the null hypothesis that metric values are equal
#             across symptomatic proportions for scenario 1.

# Post-hoc for Kruskal-Wallis
# Dunn's test performs pairwise comparisons between groups IF the Kruskal-Wallis
# test indicates a significant difference (which it doesn't in this case)
cat("\n── Dunn's Post-hoc: Ratio ──\n")

unique(R_v_I26_scenario_1$symptomatic_proportion) #"base" "base_minus_20pp" "base_plus_20pp"  "base_plus_40pp" "all_100"
unique(R_v_I26_scenario_1$symp_prop) #[1]   0 -20  20  40 100

dunn.test(R_v_I26_scenario_1$td_diff, R_v_I26_scenario_1$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 91.9776, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   0.240878
#         │     1.0000 
#         │
#      20 │   4.121219   3.884352
#         │     0.0002*    0.0005*
#         │
#      40 │   5.197029   4.961589   1.082562
#         │     0.0000*    0.0000*    1.0000 
#         │
#     100 │   7.951675   7.719159   3.841708   2.755823
#         │     0.0000*    0.0000*    0.0006*    0.0293 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)


dunn.test(R_v_I26_scenario_1$td_ratio, R_v_I26_scenario_1$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 34.0092, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │  -0.938706
#         │     1.0000 
#         │
#      20 │   1.842095   2.783120
#         │     0.3273     0.0269 
#         │
#      40 │   1.974909   2.914642   0.135658
#         │     0.2414     0.0178*    1.0000 
#         │
#     100 │   4.410352   5.351800   2.573714   2.434498
#         │     0.0001*    0.0000*    0.0503     0.0746 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)


dunn.test(R_v_I26_scenario_1$log_td_ratio, R_v_I26_scenario_1$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 34.0092, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │  -0.938706
#         │     1.0000 
#         │
#      20 │   1.842095   2.783120
#         │     0.3273     0.0269 
#         │
#      40 │   1.974909   2.914642   0.135658
#         │     0.2414     0.0178*    1.0000 
#         │
#     100 │   4.410352   5.351800   2.573714   2.434498
#         │     0.0001*    0.0000*    0.0503     0.0746 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)


# RCGP vs ICU (26 sites) for scenario 13 (winter, 500 samples per week)

cat("\n── Kruskal-Wallis: Difference ──\n")
kruskal.test(td_diff ~ symptomatic_proportion, data = R_v_I26_scenario_13)
#Kruskal-Wallis rank sum test
#data:  td_diff by symptomatic_proportion
#Kruskal-Wallis chi-squared = 136.88, df = 4, p-value < 2.2e-16

cat("\n── Kruskal-Wallis: Ratio (raw) ──\n")
kruskal.test(td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_13)
#Kruskal-Wallis rank sum test
#data:  td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 52.063, df = 4, p-value = 1.338e-10

cat("\n── Kruskal-Wallis: log(Ratio) ──\n")
kruskal.test(log_td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_13)
#Kruskal-Wallis rank sum test
#data:  log_td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 52.063, df = 4, p-value = 1.338e-10

# Post-hoc for Kruskal-Wallis
cat("\n── Dunn's Post-hoc: Ratio ──\n")
dunn.test(R_v_I26_scenario_13$td_diff, R_v_I26_scenario_13$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 136.8812, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   3.192418
#         │     0.0071*
#         │
#      20 │   5.527929   2.335511
#         │     0.0000*    0.0976 
#         │
#      40 │   7.847444   4.656623   2.322281
#         │     0.0000*    0.0000*    0.1011 
#         │
#     100 │   10.72977   7.544584   5.214357   2.894726
#         │     0.0000*    0.0000*    0.0000*    0.0190*
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)

dunn.test(R_v_I26_scenario_13$td_ratio, R_v_I26_scenario_13$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 52.0635, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   1.745780
#         │     0.4042 
#         │
#      20 │   2.860100   1.114320
#         │     0.0212*    1.0000 
#         │
#      40 │   4.322560   2.577654   1.463892
#         │     0.0001*    0.0497     0.7161 
#         │
#     100 │   6.711721   4.969891   3.858092   2.395586
#         │     0.0000*    0.0000*    0.0006*    0.0830 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)

dunn.test(R_v_I26_scenario_13$log_td_ratio, R_v_I26_scenario_13$symp_prop, method = "bonferroni")
#Kruskal-Wallis chi-squared = 52.0635, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   1.745780
#         │     0.4042 
#         │
#      20 │   2.860100   1.114320
#         │     0.0212*    1.0000 
#         │
#      40 │   4.322560   2.577654   1.463892
#         │     0.0001*    0.0497     0.7161 
#         │
#     100 │   6.711721   4.969891   3.858092   2.395586
#         │     0.0000*    0.0000*    0.0006*    0.0830 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)

# RCGP vs ICU (26 sites) for scenario 25 (winter, 1000 samples per week)

cat("\n── Kruskal-Wallis: Difference ──\n")
kruskal.test(td_diff ~ symptomatic_proportion, data = R_v_I26_scenario_25) 
#Kruskal-Wallis rank sum test
#data:  td_diff by symptomatic_proportion
#Kruskal-Wallis chi-squared = 71.004, df = 4, p-value = 1.393e-14

cat("\n── Kruskal-Wallis: Ratio (raw) ──\n")
kruskal.test(td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_25)
#Kruskal-Wallis rank sum test
#data:  td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 29.288, df = 4, p-value = 6.833e-06

cat("\n── Kruskal-Wallis: log(Ratio) ──\n")
kruskal.test(log_td_ratio ~ symptomatic_proportion, data = R_v_I26_scenario_25)
#Kruskal-Wallis rank sum test
#data:  log_td_ratio by symptomatic_proportion
#Kruskal-Wallis chi-squared = 29.288, df = 4, p-value = 6.833e-06

# Post-hoc for Kruskal-Wallis
cat("\n── Dunn's Post-hoc: Difference ──\n")
dunn.test(R_v_I26_scenario_25$td_diff, R_v_I26_scenario_25$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 71.0042, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   1.406619
#         │     0.7977 
#         │
#      20 │   4.335011   2.928392
#         │     0.0001*    0.0170*
#         │
#      40 │   5.711660   4.305393   1.377733
#         │     0.0000*    0.0001*    0.8414 
#         │
#     100 │   7.178769   5.773914   2.849193   1.472475
#         │     0.0000*    0.0000*    0.0219*    0.7045 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)

cat("\n── Dunn's Post-hoc: Ratio ──\n")
dunn.test(R_v_I26_scenario_25$td_ratio, R_v_I26_scenario_25$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 29.2876, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   0.386593
#         │     1.0000 
#         │
#      20 │   2.402203   2.015610
#         │     0.0815     0.2192 
#         │
#      40 │   3.319582   2.933086   0.917980
#         │     0.0045*    0.0168*    1.0000 
#         │
#     100 │   4.492130   4.106022   2.092939   1.175587
#         │     0.0000*    0.0002*    0.1818     1.0000 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)

cat("\n── Dunn's Post-hoc: log(Ratio) ──\n")
dunn.test(R_v_I26_scenario_25$log_td_ratio, R_v_I26_scenario_25$symp_prop, method = "bonferroni")
#Kruskal-Wallis rank sum test
#data: x and group
#Kruskal-Wallis chi-squared = 29.2876, df = 4, p-value = 0
#Dunn's Pairwise Comparison of x by group                    
#                                  (Bonferroni)                                  
#Col Mean-│
#Row Mean │        -20          0         20         40
#─────────┼────────────────────────────────────────────
#       0 │   0.386593
#         │     1.0000 
#         │
#      20 │   2.402203   2.015610
#         │     0.0815     0.2192 
#         │
#      40 │   3.319582   2.933086   0.917980
#         │     0.0045*    0.0168*    1.0000 
#         │
#     100 │   4.492130   4.106022   2.092939   1.175587
#         │     0.0000*    0.0002*    0.1818     1.0000 
#FWER = 0.05
#Reject Ho if adjusted p ≤ FWER/2, where (unadjusted) p = Pr(Z ≥ |z|)


# ════════════════════════════════════════════════════════════════════════
# STEP 6: MANOVA — TEST BOTH SIMULTANEOUSLY
# ════════════════════════════════════════════════════════════════════════
# If p-value is below significance level then cannot reject null hypothesis that
# the means are the same
cat("\n── MANOVA: Difference + Log(Ratio) jointly ──\n")
manova_model <- manova(cbind(td_diff, log_td_ratio) ~ symptomatic_proportion, data = R_v_I26_scenario_1)
summary(manova_model, test = "Pillai")   # Pillai's trace (most robust)
#                         Df   Pillai approx F num Df den Df    Pr(>F)    
#symptomatic_proportion    4 0.087162   56.617      8   9940 < 2.2e-16 ***
#  Residuals              4970                                              
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(manova_model, test = "Wilks")    # Wilks' lambda (most common)
#                         Df   Wilks approx F num Df den Df    Pr(>F)    
#symptomatic_proportion    4 0.91303   57.816      8   9938 < 2.2e-16 ***
#  Residuals              4970                                             
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Follow-up univariate ANOVAs from MANOVA
summary.aov(manova_model)
#Response td_diff :
#                         Df Sum Sq Mean Sq F value    Pr(>F)    
#symptomatic_proportion    4   7470 1867.52  22.785 < 2.2e-16 ***
#  Residuals              4970 407344   81.96                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Response log_td_ratio :
#                         Df  Sum Sq Mean Sq F value    Pr(>F)    
#symptomatic_proportion    4   1.606 0.40148  7.3223 7.088e-06 ***
#  Residuals              4970 272.501 0.05483                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#2 observations deleted due to missingness


# ════════════════════════════════════════════════════════════════════════
# STEP 7: VISUALIZE
# ════════════════════════════════════════════════════════════════════════

# Colour blind friendly palette
cbPalette <- c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9")

# Function to change symptomatic proportion labels so appear in numeric order in plots
library(dplyr)
add_plot_codes <- function(df){
  df <- df %>%
    mutate(symp_plot_code = case_when(
      symptomatic_proportion == "base_minus_20pp" ~ "A",
      symptomatic_proportion == "base" ~ "B",
      symptomatic_proportion == "base_plus_20pp" ~ "C",
      symptomatic_proportion == "base_plus_40pp" ~ "D",
      symptomatic_proportion == "all_100" ~ "E",
    ))
  return(df)
}

# Function to plot
plot_diff_ratio_boxplot <- function(df){
  ggplot(df, aes(x = symp_plot_code, y = value, fill = symp_plot_code)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = cbPalette) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~metric, scales = "free_y",
               labeller = as_labeller(c(td_diff  = "Difference (TD_RCGP - TD_ICU_26_sites)",# (D = X - Y)",
                                        #log_td_ratio   = "Log(TD_RCGP / TD_ICU_26_sites)"))) + # (log(X/Y))"))) +
                                        td_ratio   = "TD_RCGP / TD_ICU_26_sites"))) + # (log(X/Y))"))) +
    labs(title = "Effect of symptomatic proportion on Difference and Ratio of TDs",
         x = "Symptomatic proportion", y = "") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")
}

plot_boxplots <- function(df, yval, ylimits){
  ggplot(df, aes(x = symp_plot_code, y = yval, fill = symp_plot_code)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = cbPalette) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    coord_cartesian(ylim = ylimits) +
    facet_wrap(~scenario_number, scales = "free_y",
               labeller = as_labeller(c("1"="300 samples per week","13"="500","25"="1000"))) + 
    labs(title = "Effect of symptomatic proportion on Difference and Ratio of TDs",
         x = "Symptomatic proportion", y = "") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")
  
}


### Scenario 1 (300 samples per week, winter)
data_long_RvI26_scen1 <- pivot_longer(R_v_I26_scenario_1, cols = c(td_diff, td_ratio),#log_td_ratio),
                                    names_to = "metric", values_to = "value")
data_long_RvI26_scen1 <- add_plot_codes( df = data_long_RvI26_scen1)
unique(data_long_RvI26_scen1$symptomatic_proportion) #[1] "base" "base_minus_20pp" "base_plus_20pp"  "base_plus_40pp"  "all_100"
p_RvI26_scen1 <- plot_diff_ratio_boxplot(df=data_long_RvI26_scen1)
ggsave( paste0(base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_diff_ratio_scenario_1.jpeg"), plot = p_RvI26_scen1)

# Scenario 13 (500 samples per week, winter)
data_long_RvI26_scen13 <- pivot_longer(R_v_I26_scenario_13, cols = c(td_diff, log_td_ratio),
                                     names_to = "metric", values_to = "value")
data_long_RvI26_scen13 <- add_plot_codes( df = data_long_RvI26_scen13 )
p_RvI26_scen13 <- plot_diff_ratio_boxplot( df = data_long_RvI26_scen13 )
ggsave( paste0( base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_diff_ratio_scenario_13.jpeg"), plot = p_RvI26_scen13)

# Scenario 25 (1000 samples per week, winter)
data_long_RvI26_scen25 <- pivot_longer(R_v_I26_scenario_25, cols = c(td_diff, log_td_ratio),
                                     names_to = "metric", values_to = "value")
data_long_RvI26_scen25 <- add_plot_codes( df = data_long_RvI26_scen25 )
p_RvI26_scen25 <- plot_diff_ratio_boxplot( df = data_long_RvI26_scen25 )
ggsave( paste0( base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_diff_ratio_scenario_25.jpeg"), plot = p_RvI26_scen25)


## Manipulate df to try to plot all scenarios together
rcgp_v_icu26 <- subset( td_symp_diff_ratio_df, sampling_methods =="RCGP_v_ICU_sites_26" & !is.infinite(td_diff))
rcgp_v_icu26 <- rcgp_v_icu26[,c("symptomatic_proportion","scenario_number","sim_rep","td_diff","td_ratio")]
rcgp_v_icu26 <- add_plot_codes(df=rcgp_v_icu26)

# Check min and max values for plots
rcgp_v_icu26_filtered_td_diff <- rcgp_v_icu26 %>% filter(!is.nan(td_diff))
min( rcgp_v_icu26_filtered_td_diff$td_diff ) #[1] -40.40725
max( rcgp_v_icu26_filtered_td_diff$td_diff ) #[1] 56.24702
rcgp_v_icu26_filtered_td_ratio <- rcgp_v_icu26 %>% filter(!is.nan(td_ratio))
min( rcgp_v_icu26_filtered_td_ratio$td_ratio ) #[1] 0.1529307
max( rcgp_v_icu26_filtered_td_ratio$td_ratio ) #[1] 7.935463


# Difference in TD
p_RvI26_diff_all_scenarios <- 
  ggplot(rcgp_v_icu26, aes(x = symp_plot_code, y = td_diff, fill = symp_plot_code)) +
    geom_jitter(width = 0.15, alpha = 0.2) +  
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = cbPalette) +
    coord_cartesian(ylim = c(-40,+60)) + # See min and max values above
    facet_wrap(~scenario_number, scales = "free_y",
               labeller = as_labeller(c("1"="300 samples per week","13"="500","25"="1000"))) + 
    labs(title = "Effect of symptomatic proportion on Difference and Ratio of TDs",
         x = "Symptomatic proportion", y = "") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")
p_RvI26_diff_all_scenarios
ggsave( paste0( base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_diff_all_scenarios.jpeg"), plot = p_RvI26_diff_all_scenarios)
# Also save in Supplementary Material Figures folder
figs_folder = dirname(dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))))
ggsave( paste0( figs_folder, "/results/figures_supplementary/raw/RCGP_v_ICU26_diff_all_scenarios.jpeg"), plot = p_RvI26_diff_all_scenarios)

# Ratio of TD
p_RvI26_ratio_all_scenarios <- 
ggplot(rcgp_v_icu26, aes(x = symp_plot_code, y = td_ratio, fill = symp_plot_code)) +
  geom_jitter(width = 0.15, alpha = 0.2) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = cbPalette) +
  coord_cartesian(ylim = c(0,8)) + # See min and max values above
  facet_wrap(~scenario_number, scales = "free_y",
             labeller = as_labeller(c("1"="300 samples per week","13"="500","25"="1000"))) + 
  labs(title = "Effect of symptomatic proportion on Difference and Ratio of TDs",
       x = "Symptomatic proportion", y = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

p_RvI26_ratio_all_scenarios

ggsave( paste0( base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_ratio_all_scenarios.jpeg"), plot = p_RvI26_ratio_all_scenarios)
# Also save in Supplementary Material Figures folder
ggsave( paste0( figs_folder, "/results/figures_supplementary/raw/RCGP_v_ICU26_ratio_all_scenarios.jpeg"), plot = p_RvI26_ratio_all_scenarios)

## Manipulate df to try to plot all scenarios together
rcgp_v_icu26_log_td_ratio <- subset( td_symp_diff_ratio_df, sampling_methods =="RCGP_v_ICU_sites_26" & !is.infinite(td_diff))
rcgp_v_icu26_log_td_ratio <- rcgp_v_icu26_log_td_ratio[,c("symptomatic_proportion","scenario_number","sim_rep","log_td_ratio")]
rcgp_v_icu26_log_td_ratio <- add_plot_codes(df=rcgp_v_icu26_log_td_ratio)

# Check min and max values for plots
rcgp_v_icu26_log_td_ratio_filtered_td_diff <- rcgp_v_icu26_log_td_ratio %>% filter(!is.nan(log_td_ratio))
min( rcgp_v_icu26_log_td_ratio_filtered_td_diff$log_td_ratio ) #[1] -1.87777
max( rcgp_v_icu26_log_td_ratio_filtered_td_diff$log_td_ratio ) #[1] 2.071342

# Log (ratio of TDs)
p_RvI26_log_ratio_all_scenarios <- 
  ggplot(rcgp_v_icu26_log_td_ratio, aes(x = symp_plot_code, y = log_td_ratio, fill = symp_plot_code)) +
  geom_jitter(width = 0.15, alpha = 0.2) +  
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = cbPalette) +
  coord_cartesian(ylim = c(-2,+2.5)) + # See min and max values above
  facet_wrap(~scenario_number, scales = "free_y",
             labeller = as_labeller(c("1"="300 samples per week","13"="500","25"="1000"))) + 
  labs(title = "Effect of symptomatic proportion on Difference and Log of Ratio of TDs",
       x = "Symptomatic proportion", y = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
p_RvI26_log_ratio_all_scenarios
ggsave( paste0( base_folder, "/covid_like_varied_sympt_prop/RCGP_v_ICU26_log_ratio_all_scenarios.jpeg"), plot = p_RvI26_log_ratio_all_scenarios)
# Also save in Supplementary Material Figures folder
ggsave( paste0( figs_folder, "/results/figures_supplementary/raw/RCGP_v_ICU26_log_ratio_all_scenarios.jpeg"), plot = p_RvI26_log_ratio_all_scenarios)