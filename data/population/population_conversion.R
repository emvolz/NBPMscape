# Identify population by ITL2 geo area

library(readxl)

# Read population data
#'mye22final.xslx' is mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. Source: ONS
#'https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland
#[Accessed 6 November 2024]
excel_sheets("~/mSCAPE/1_epidemic_modelling/population/mye22final.xlsx")
mye22final_pop <- read_excel("~/mSCAPE/1_epidemic_modelling/population/mye22final.xlsx", sheet = "MYEB1")

unique(mye22final_pop$geography)
# [1] "Country"                   "Region"                    "Unitary Authority"         "Metropolitan County"      
# [5] "Metropolitan District"     "County"                    "Non-metropolitan District" "London Borough"           
# [9] "Council Area"              "Local Government District"


# Read files for converting between geographies
LAD2024_geo_conv <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD_(December_2024)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2025)_Lookup_in_the_UK.csv")
excel_sheets("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD20_LAU121_ITL321_ITL221_ITL121_UK_LU_v2.xlsx")
LAD20_LAU121_ITL21_UK_LU_v2 <- read_excel("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD20_LAU121_ITL321_ITL221_ITL121_UK_LU_v2.xlsx"
                                          , sheet = "LAD20_LAU121_ITL21_UK_LU v2")

# Try to match the geo areas 

mye22final_pop_local_gov_district <- unique( as.data.frame( subset( mye22final_pop, mye22final_pop$geography == "Local Government District" ) )$code )
mye22final_pop_unitary_authority <- unique( as.data.frame( subset( mye22final_pop, mye22final_pop$geography == "Unitary Authority" ) )$code )
mye22final_pop_council_area <- unique( as.data.frame( subset( mye22final_pop, mye22final_pop$geography == "Council Area" ) )$code )

sum( mye22final_pop_local_gov_district %in% unique( LAD2024_geo_conv$LAD24CD ) ) # [1] 11

table( mye22final_pop_unitary_authority %in% unique( LAD2024_geo_conv$LAD24CD ) ) # [1] TRUE 85

length(mye22final_pop_unitary_authority) # [1] 85
length( unique( LAD2024_geo_conv$LAD24CD ) ) # [1] 361

table( unique( mye22final_pop$code ) %in% unique( LAD2024_geo_conv$LAD24CD ) ) # [1] FALSE 43 TRUE 361
# Looks like all LAD24CD are in the population data

'%ni%' <- Negate('%in%')

# Filter population data for codes that match an LAD24CD
keep_index <- which( mye22final_pop$code %in% unique( LAD2024_geo_conv$LAD24CD ) )
mye22final_pop_filtered <- mye22final_pop[ keep_index , ]
# Does the population when filtered by the LAD24CD match the sum for country level
# (i.e. make sure we're not losing any population by filtering)
sum( mye22final_pop_filtered$population_2022 ) # [1] 67596281
sum( subset( mye22final_pop, mye22final_pop$name == "UNITED KINGDOM" )$population_2022 ) # [1] 67596281
# Numbers match so filtering by LAD24CD retains all the population

# Add ITL2 code and name to population data
mye22final_pop_filtered_w_ITL2 <- merge( x = mye22final_pop_filtered
                                         , y = LAD2024_geo_conv[ , c("ITL225CD","ITL225NM","LAD24CD","LAD24NM") ]
                                        , by.x = "code", by.y = "LAD24CD"
                                        , all.x = TRUE, all.y = FALSE
                                         )
# Check
sum( mye22final_pop_filtered_w_ITL2$population_2022 ) # [1] 67996751 !!! POPULATION HAS INCREASED !!!
# Try different method
mye22final_pop_filtered_w_ITL2 <- merge( x = mye22final_pop_filtered
                                         , y = LAD2024_geo_conv[ , c("ITL225CD","ITL225NM","LAD24CD","LAD24NM") ]
                                         , by.x = "code", by.y = "LAD24CD"
                                         , all.x = TRUE, all.y = FALSE
)
# Check
sum( mye22final_pop_filtered_w_ITL2$population_2022 ) # [1] 67729771 !!! POPULATION HAS INCREASED !!!
# Try a different method
mye22final_pop_filtered_w_ITL2 <- mye22final_pop_filtered
mye22final_pop_filtered_w_ITL2[,c("ITL225CD","ITL225NM")] <- NA
for (i in 1:nrow( mye22final_pop_filtered_w_ITL2 ) ){
  code_i <- mye22final_pop_filtered_w_ITL2$code[ i ]
  mye22final_pop_filtered_w_ITL2$ITL225CD[ i ] <- subset( LAD2024_geo_conv, LAD2024_geo_conv$LAD24CD == code_i )$ITL225CD
  mye22final_pop_filtered_w_ITL2$ITL225NM[ i ] <- subset( LAD2024_geo_conv, LAD2024_geo_conv$LAD24CD == code_i )$ITL225NM
}
# Check
sum( mye22final_pop_filtered_w_ITL2$population_2022 ) # [1] 67596281 # Population total matches previous totals. Number of rows also match

mye22final_pop_filtered_ITL2_age_sex <- mye22final_pop_filtered_w_ITL2[, c("ITL225CD","ITL225NM","sex","age","population_2022")]

library(dplyr)

mye22final_pop_filtered_ITL2_age_sex_collapsed <- mye22final_pop_filtered_ITL2_age_sex %>%
                                                    group_by(ITL225CD, ITL225NM, sex, age) %>%
                                                    summarize(total_population_2022 = sum(population_2022)) %>%
                                                    ungroup()
# Check
sum( mye22final_pop_filtered_ITL2_age_sex_collapsed$total_population_2022 ) # [1] 67596281 # Total population unchanged

# Version without disaggregation by sex and age
mye22final_pop_filtered_ITL2_collapsed <- mye22final_pop_filtered_ITL2_age_sex_collapsed[,c("ITL225CD","ITL225NM","total_population_2022")]
mye22final_pop_filtered_ITL2_collapsed <- mye22final_pop_filtered_ITL2_collapsed %>%
                                            group_by(ITL225CD, ITL225NM) %>%
                                            summarize(total_population_2022 = sum(total_population_2022)) %>%
                                            ungroup()
# Check
sum( mye22final_pop_filtered_ITL2_collapsed$total_population_2022 ) # [1] 67596281 # Total population unchanged

# Save data frames to list
ITL2_pop_MYE_2022 <- list(  as.data.frame( mye22final_pop_filtered_ITL2_collapsed )
                          , as.data.frame( mye22final_pop_filtered_ITL2_age_sex_collapsed )
                          )
# Change names of items in list
names( ITL2_pop_MYE_2022 )[1] <- "population"
names( ITL2_pop_MYE_2022 )[2] <- "population_by_age_sex"

# Save file
saveRDS( ITL2_pop_MYE_2022 , "~/mSCAPE/1_epidemic_modelling/population/ITL2_pop_MYE_2022.rds" )
# Check
test <- readRDS("~/mSCAPE/1_epidemic_modelling/population/ITL2_pop_MYE_2022.rds" )
View(test[1]$population)


# Plot map

# Load required libraries
library(ggplot2)
library(sf)
library(dplyr)

# Load the shapefile for ITL2 regions (replace 'path_to_shapefile' with the actual file path)
# https://ckan.publishing.service.gov.uk/dataset/international-territorial-level-2-january-2025-boundaries-uk-bgc/resource/78abd24e-5f15-4e04-a118-ddc649e9c769
itl2_shapefile <- st_read("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2 shape files/ITL2_JAN_2025_UK_BFE.shp")

# Example data frame with ITL2 region codes and values to color by
region_data <- data.frame( ITL225CD = ITL2_pop_MYE_2022$population$ITL225CD # ITL2 codes
                        , pop_mye_2022 = ITL2_pop_MYE_2022$population$total_population_2022 # ONS Census population estimates mid-year 2022
)

# Join the shapefile data with the region data
itl2_map_data <- itl2_shapefile %>%
  left_join( region_data, by = c("ITL225CD" = "ITL225CD")) #c("ITL2_Code" = "ITL2_Code"))

# Plot the map using ggplot2
ggplot(data = itl2_map_data) +
  geom_sf(aes(fill = pop_mye_2022), color = "black") + # Fill regions by 'pax_2024_per_day_integer' column
  scale_fill_gradient(low = "lightgreen", high = "darkgreen", na.value = "white") + # Adjust color gradient
  theme_minimal() +
  labs(
    title = "UK Map with ITL2 Regions",
    fill = "ONS population estimte mid-year 2022"
  )
