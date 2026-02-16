# Script to obtain data on international travel to the UK
# specifically age profile and regions visited/return to

library(readxl)

# ONS travelpac data for 2009-2019 (there is an issue with later data)
# Source: https://www.ons.gov.uk/peoplepopulationandcommunity/leisureandtourism/datasets/travelpac
# Downloaded 2 May 2025
travelpac <- read_excel("~/mSCAPE/1_epidemic_modelling/international_travel/travelpac2009_2019/TravelPac 2009-2019 Labelled.xlsx", sheet = "Labelled 2009-2019")
table( travelpac$ukos )

# Aggregate number of international traveller visits by residency, and age group
travelpac_2019 <- subset( travelpac , travelpac$Year == 2019 )
travelpac_2019_UK <- subset( travelpac_2019 , travelpac_2019$ukos == "UK residents" )
travelpac_2019_OS <- subset( travelpac_2019 , travelpac_2019$ukos == "Overseas residents" )
sum( travelpac_2019$visits , na.rm = TRUE ) #[1] 133,943,220
sum( travelpac_2019_UK$visits , na.rm = TRUE ) #[1] 93,085,772
sum( travelpac_2019_OS$visits , na.rm = TRUE ) #[1] 40,857,448

travelpac_2019_UK_aggregated <- aggregate(visits ~ Age, data=travelpac_2019_UK, FUN=sum)
sum( travelpac_2019_UK_aggregated$visits ) #[1] 93,085,772

travelpac_2019_OS_aggregated <- aggregate(visits ~ Age, data=travelpac_2019_OS, FUN=sum)
sum( travelpac_2019_OS_aggregated$visits ) #[1] 40,857,448

travelpac_2019_UK_aggregated <- cbind( replicate(8,"UK residents") , travelpac_2019_UK_aggregated )
colnames( travelpac_2019_UK_aggregated )[1] <- "ukos"

travelpac_2019_OS_aggregated <-cbind( replicate(8,"Overseas residents") , travelpac_2019_OS_aggregated )
colnames( travelpac_2019_OS_aggregated )[1] <- "ukos"

travelpac_2019_aggregated <- rbind( travelpac_2019_UK_aggregated, travelpac_2019_OS_aggregated )

# Allocate visits with unknown age group to other age groups based on known age group weights
uk_weights <- travelpac_2019_UK_aggregated[ , "visits"] / sum( travelpac_2019_UK_aggregated[ , "visits"] ) 
os_weights <- travelpac_2019_OS_aggregated[ , "visits"] / sum( travelpac_2019_OS_aggregated[ , "visits"] ) 

uk_unknown_age_visits <- subset( travelpac_2019_UK_aggregated , travelpac_2019_UK_aggregated$Age == "D/K" )["visits"]
os_unknown_age_visits <- subset( travelpac_2019_OS_aggregated , travelpac_2019_OS_aggregated$Age == "D/K" )["visits"]

travelpac_2019_UK_aggregated_inc_unknown <- travelpac_2019_UK_aggregated[1:7,]
travelpac_2019_UK_aggregated_inc_unknown[,"visits"] <- travelpac_2019_UK_aggregated_inc_unknown[,"visits"] + 
                                                              ( mapply( '*', uk_unknown_age_visits, uk_weights[1:7] ) )

travelpac_2019_OS_aggregated_inc_unknown <- travelpac_2019_OS_aggregated[1:7,]
travelpac_2019_OS_aggregated_inc_unknown[,"visits"] <- travelpac_2019_OS_aggregated_inc_unknown[,"visits"] + 
                                                                ( mapply( '*', os_unknown_age_visits, os_weights[1:7] ) )

# Combine UK and Overseas resident visits
travelpac_2019_adjusted <- travelpac_2019_OS_aggregated_inc_unknown
travelpac_2019_adjusted[,"ukos"] <- NULL
travelpac_2019_adjusted[,"visits"] <- travelpac_2019_UK_aggregated_inc_unknown[,"visits"] + travelpac_2019_OS_aggregated_inc_unknown[,"visits"]
#Check
sum( travelpac_2019_adjusted[,"visits"] ) == sum( as.data.frame( travelpac_2019 )[,"visits"] , na.rm = T )
#**FALSE**
# 133,361,808 != 133,943,220
sum( travelpac_2019_adjusted[,"visits"] ) == sum( travelpac_2019_UK_aggregated[,"visits"] ) + sum( travelpac_2019_OS_aggregated[,"visits"] )
#**FALSE**
# 133,361,808 != (133,943,220 = 93,085,772 + 40,857,448) 

travelpac_2019_age_group_weights <- travelpac_2019_adjusted
travelpac_2019_age_group_weights[,"percent"] <- travelpac_2019_adjusted[, "visits"] / sum( travelpac_2019_adjusted[, "visits"] ) 

# Save output
saveRDS( travelpac_2019_age_group_weights, "~/mSCAPE/1_epidemic_modelling/international_travel/travelpac2009_2019/travelpac_2019_age_group_weights.rds" )


# Allocate age group weights to single year age using ONS UK population data

# Read in population
# 'mye22final.xslx' is mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. Source: ONS
# [Accessed 6 November 2024]
# Available at #'https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
ons_pop <- read_excel("~/mSCAPE/1_epidemic_modelling/population/mye22final.xlsx", sheet = "MYEB1")
unique(ons_pop$name)
unique(ons_pop$sex)
# Filter for UK geo level
ons_pop_uk <- as.data.frame( subset(ons_pop, ons_pop$name == "UNITED KINGDOM") )
# Combine male and female (because not currently differentiating in NBPMscape model)
ons_pop_uk_all_sex <- aggregate(population_2022 ~ age, data = ons_pop_uk, sum)
#Check
sum(ons_pop_uk_all_sex$population_2022) #[1] 67596281
sum(ons_pop_uk$population_2022)         #[1] 67596281

# Add percent column
ons_pop_uk_all_sex$proportion <- ons_pop_uk_all_sex$population_2022 / sum( ons_pop_uk_all_sex$population_2022 )

# Allocate age group weights to single year ages
parse_range <- function(str) {
  # Split string range into start and end age of age group
  #parts <- as.integer( strsplit(str, "-")[[1]] )
  # If the age group name doesn't follow the standard format (e.g. 65 & over) then adjust
  # Assume it is the final age group and use the maximum age in UK ONS population data (90 years)
  if ( is.na( as.integer( strsplit(str, "-")[[1]][1] ) ) ){
    parts <- as.integer( c( as.integer( strsplit(str, " ")[[1]][1] ) , 90 ) )
  } else {
    parts <- as.integer( strsplit(str, "-")[[1]] )
  }
  return( seq( parts[1], parts[2] ) )
}

# Cycle through age groups and convert to single year age proportions
age_group_single_years <- as.data.frame( matrix( data = NA, nrow = 1, ncol = 4 
                                                 , byrow = FALSE, dimnames = NULL ) )
colnames( age_group_single_years ) <- c( colnames( ons_pop_uk_all_sex ) , "intl_travel_prop" )
for (i in 1: nrow( travelpac_2019_age_group_weights ) ){
  n_digits <- nchar(gsub("[^0-9]", "", travelpac_2019_age_group_weights$percent[ i ]))
  temp_single_years <- parse_range( travelpac_2019_age_group_weights$Age[ i ] )
  age_group_weight <- travelpac_2019_age_group_weights$percent[ i ]
  pop_single_age_dist <- subset( ons_pop_uk_all_sex 
                                 , ons_pop_uk_all_sex$age %in% temp_single_years )
  pop_single_age_dist$intl_travel_prop <- age_group_weight *
                                          ( pop_single_age_dist$proportion /
                                             sum( pop_single_age_dist$proportion ) )
  
  # Check haven't lost any of the weight
  if( round( sum( pop_single_age_dist$intl_travel_prop ),7) == round( age_group_weight, 7) ){
  } else{ print( c( "Error: weights do not match for ", travelpac_2019_age_group_weights$Age[i])) }
  
  
  age_group_single_years <- rbind( age_group_single_years, pop_single_age_dist )
}
# Remove empty first line
if( sum( is.na(age_group_single_years[1,]) ) == 4 ){
  age_group_single_years <- age_group_single_years[2:nrow(age_group_single_years),]
}
rownames( age_group_single_years ) <- NULL

# Check haven't lost any of the weight
sum( age_group_single_years$intl_travel_prop ) #[1] 1
sum( travelpac_2019_age_group_weights$percent ) #[1] 1


# Plot shows that there are clear issues with this approach on its own.
# Large step changes and also the 90+ weight allocated to single year ago = 90
plot( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop , typ  = "o" , xlab = "single year age", ylab = "proportion of population" )
lines( x = age_group_single_years$age , y = age_group_single_years$proportion , typ  = "o", col = "red"  )
legend("topright"
       , legend = c("Int'l travellers adjusted"# for UK population distribution"
                    , "UK population")
       , col = c("black","red")
       , lty = 1#c("l","l")
       )

# Try to adjust for discontinuities in adjusted data
# Identify change points in data and smooth distribution
install.packages("changepoint")
library(changepoint)
cpt <- cpt.mean( age_group_single_years$intl_travel_prop
                 , penalty = "None"
                 , method = "BinSeg" #"PELT"#
                 , Q = 5 ) # Detects mean shifts (steps)
change_points <- cpts(cpt)

# Adjust the data to remove steps:
cpt_adjusted_data <- age_group_single_years$intl_travel_prop
for (i in seq_along(change_points)) {
  if (i == 1) next
  step_size <- mean( age_group_single_years$intl_travel_prop[(change_points[i-1]+1):change_points[i]]) - 
                                           mean(age_group_single_years$intl_travel_prop[1:change_points[i-1]])
  cpt_adjusted_data[(change_points[i]+1):length(age_group_single_years$intl_travel_prop)] <- 
    cpt_adjusted_data[(change_points[i]+1):length(age_group_single_years$intl_travel_prop)] - step_size
}       

# Plot shows that it doesn't seem to work
plot( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop , typ  = "o" , xlab = "single year age", ylab = "proportion of population" )
lines( x = age_group_single_years$age , y = age_group_single_years$proportion , typ  = "o", col = "red"  )
lines( x = age_group_single_years$age , y = cpt_adjusted_data , typ  = "o", col = "blue"  )
legend("topright"
       , legend = c("Int'l travellers adjusted"# for UK population distribution"
                    , "UK population"
                    , "Int'l travellers changepoint adjusted")
       , col = c("black","red", "blue")
       , lty = 1#c("l","l")
)

# Leave with discontinuities but adjust for 90+ age group in UK ONS population data
# Assume that no international travel for individuals aged 90+
age_group_single_years[92:101,] <- data.frame(age = seq(91,100)
                                              , population_2022 = rep(NA,10)
                                              , proportion = rep(NA,10)
                                              , intl_travel_prop = rep(NA,10))
age_group_single_years$intl_travel_prop_adj <- age_group_single_years$intl_travel_prop

# Different options for redistributing the 90+ weight
# 1) Equal weight
# 2) Linear decline between 90 and 100
# 3) exponential decay between 90 and 100

# 1) Equal weight to between 90 and 100 years
age_group_single_years$intl_travel_prop_adj[91:101] <- age_group_single_years$intl_travel_prop[91] / 11
sum(age_group_single_years$intl_travel_prop_adj[91:101]) == age_group_single_years$intl_travel_prop[91]
# Plot
plot( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop , typ  = "o" , xlab = "single year age", ylab = "proportion of population" )
lines( x = age_group_single_years$age , y = age_group_single_years$proportion , typ  = "o", col = "red"  )
lines( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop_adj , typ  = "o", col = "blue"  )
legend("topright"
       , legend = c("Int'l travellers adjusted"# for UK population distribution"
                    , "UK population"
                    , "Int'l travellers adjusted for 90+ age")
       , col = c("black","red", "blue")
       , lty = 1#c("l","l")
)

# 2) Linear decline to zero at 101 years
weight_to_redistribute <- age_group_single_years$intl_travel_prop[91]
equal_increments <- weight_to_redistribute / 11
t <- rev( seq(90,101)) - 90
t_sum <- sum(t)
linear_increments <- weight_to_redistribute / t_sum
linear_increments_vec <- linear_increments * t
age_group_single_years$intl_travel_prop_adj[91:101] <- linear_increments_vec[1:11]

sum(age_group_single_years$intl_travel_prop_adj[91:101]) == age_group_single_years$intl_travel_prop[91]

# Plot
plot( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop , typ  = "o" , xlab = "single year age", ylab = "proportion of population", ylim=c(0,0.025) , col = "red", pch=16)
lines( x = age_group_single_years$age , y = age_group_single_years$proportion , typ  = "o", col = "black" , pch=16 )
lines( x = age_group_single_years$age , y = age_group_single_years$intl_travel_prop_adj , typ  = "o", col = "blue", pch = 16  )
legend("topright"
       , legend = c("UK population"
                    , "Int'l travellers by age group"
                    , "Int'l travellers gen pop adjusted"# for UK population distribution"
                    , "Int'l travellers gen pop adj (90+ age)")
       , col = c("black","green", "red", "blue")
       , lty = 1#c("l","l")
       , pch = c(16,NA,16,1)
       , cex = 0.9
       , lwd = 2
)
abline(h=0)
lines( x = seq(0,15), y = rep(travelpac_2019_age_group_weights$percent[1]/(15-0+1),(15-0+1)), col="green", lwd=2)
lines( x = seq(16,24), y = rep(travelpac_2019_age_group_weights$percent[2]/(24-16+1),(24-16+1)), col="green", lwd=2)
lines( x = seq(25,34), y = rep(travelpac_2019_age_group_weights$percent[3]/(34-25+1),(34-25+1)), col="green", lwd=2)
lines( x = seq(35,44), y = rep(travelpac_2019_age_group_weights$percent[4]/(44-35+1),(44-35+1)), col="green", lwd=2)
lines( x = seq(45,54), y = rep(travelpac_2019_age_group_weights$percent[5]/(54-45+1),(54-45+1)), col="green", lwd=2)
lines( x = seq(55,64), y = rep(travelpac_2019_age_group_weights$percent[6]/(64-55+1),(64-55+1)), col="green", lwd=2)
lines( x = seq(65,100), y = rep(travelpac_2019_age_group_weights$percent[7]/(100-65+1),(100-65+1)), col="green", lwd=2)

# Potential issues with this approach
# - Quite a lot of adjustment to avoid the step changes in the age distribution 
#   for international travel
# - assumption that distribution of international travellers is the same as general 
#   population by age group may not be correct - particularly likely at the higher 
#   age end given the breadth of the oldest category (65 & older)
# - mismatch between international traveller data (2019) and population data (2022)
# - may be an impact from excess deaths during pandemic on population distribution
#   and international travel behaviour within certain age groups

saveRDS( age_group_single_years, "~/mSCAPE/1_epidemic_modelling/international_travel/travelpac2009_2019/travelpac_2019_age_single_year_weights.rds" )
