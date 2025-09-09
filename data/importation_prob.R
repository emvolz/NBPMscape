# Compute importation probability by region based on passenger numbers
# using Civil Aviation Authority data
# March 2025

# 1 Load CAA data
# 2 Manipulate CAA data
# 3 Assign Region, UTLA and LTLA to airport

# Note that arrivals and departures are not differentiated (also probably includes stops/connections)


################################################################################
# 1  Load CAA data
################################################################################
# https://www.caa.co.uk/data-and-analysis/uk-aviation-market/airports/uk-airport-data/uk-airport-data-2024/
# Table 12.1 is the dataset used

# Read in data
# see 'CAA_data_download.R' for full script for downloading and reading in files
# but here we use part of this script

# Read in files
################################################################################

# Create lists of files to read
folder_path <- "~/mSCAPE/1_epidemic_modelling/Air traffic/Source data/Civil Aviation Authority/"

# Files for 2024
file_list_table_12_1 <- list.files( path = folder_path , pattern = "^Table_12_1.*2024.*\\.csv$" )

# Create dataframe to store and append the data read from files
caa_table_12_1_df <- as.data.frame( data.table::fread( paste0( folder_path, file_list_table_12_1[ 1 ] ) ) )

# Read in table 12.1 files and append data to dataframe
for( f in 2:length( file_list_table_12_1 ) ){
  
  df_name <- stringr::str_split( file_list_table_12_1[ f ], pattern = ".csv" )[[1]][1]
  
  df_temp <- data.table::fread( paste0( folder_path, file_list_table_12_1[ f ] )
                                #, sep = "\c"
                                #, drop = c("Pax Share","Fare","Est.","Revenue","Yield","RPM")
                                #, blank.lines.skip = TRUE
  )
  df_temp <- as.data.frame( df_temp )
  
  # Check columns names match - necessary for rbind to work
  check_df_colnames <- mean( colnames( caa_table_12_1_df ) == colnames( df_temp ) )
  print( paste( file_list_table_12_1[ f ] , check_df_colnames ) )
  
  # Append to combined dataframe
  caa_table_12_1_df <- rbind( caa_table_12_1_df , df_temp )
  rm( df_temp )
}

# Check
HEATHROW    <- sum( subset( caa_table_12_1_df, UK_airport == "HEATHROW" )$total_pax_this_period ) / 2
#[1] 39414917
MANCHESTER <- sum( subset( caa_table_12_1_df, UK_airport == "MANCHESTER" )$total_pax_this_period ) / 2
#[1] 14427330
HEATHROW/MANCHESTER
#[1] 2.731962


################################################################################
# 2 Manipulate CAA data
################################################################################
# Remove columns that are not of interest
caa_table_12_1_df[,c("total_pax_scheduled_this_period"
                    ,"total_pax_charter_this_period"
                    ,"total_pax_last_period"
                    ,"total_pax_scheduled_last_period"
                    ,"total_pax_charter_last_period"
                    ,"total_pax_percent")] <- NULL

# Create 3D array (UK airport x non-UK airport (or country or larger foreign region) x month in 2024)

uk_airports       <- sort( unique( caa_table_12_1_df$UK_airport ) )
foreign_regions   <- sort( unique( caa_table_12_1_df$foreign_region ) )
foreign_countries <- sort( unique( caa_table_12_1_df$foreign_country ) )
foreign_airports  <- sort( unique( caa_table_12_1_df$foreign_airport ) )

# Create the 3D arrays
airport_airport_pax_array   <- array(0, dim = c( length( uk_airports ), length( foreign_airports  ), 12 ) )
airport_countries_pax_array <- array(0, dim = c( length( uk_airports ), length( foreign_countries ), 12 ) )
airport_regions_pax_array   <- array(0, dim = c( length( uk_airports ), length( foreign_regions   ), 12 ) )

# Fill arrays with summed values
# airport_airport_array
for (i in 1:nrow( caa_table_12_1_df ) ) {
  row_index <- match( caa_table_12_1_df[i, "UK_airport"], uk_airports )
  col_index <- match( caa_table_12_1_df[i, "foreign_airport"], foreign_airports )
  depth_index <- as.numeric( stringr::str_sub( ( caa_table_12_1_df[ i, "this_period"] ), -2) )
  airport_airport_pax_array[row_index, col_index, depth_index] <- airport_airport_pax_array[row_index, col_index, depth_index] +
                                                          caa_table_12_1_df[i, "total_pax_this_period"]
}
dimnames( airport_airport_pax_array ) <- list( uk_airports, foreign_airports, seq(1,12,1))

# Check
View( airport_airport_pax_array[ "GATWICK",,])
GATWICK <- sum( airport_airport_pax_array[ "GATWICK",,]) / 2    #[1] 20139357
MANCHESTER <- sum( airport_airport_pax_array[ "MANCHESTER",,]) / 2 #[1] 14427330
GATWICK/MANCHESTER #[1] 1.395917

# airport_countries_array
for (i in 1:nrow( caa_table_12_1_df ) ) {
  row_index <- match( caa_table_12_1_df[i, "UK_airport"], uk_airports )
  col_index <- match( caa_table_12_1_df[i, "foreign_country"], foreign_countries )
  depth_index <- as.numeric( stringr::str_sub( ( caa_table_12_1_df[ i, "this_period"] ), -2) )
  airport_countries_pax_array[row_index, col_index, depth_index] <- airport_countries_pax_array[row_index, col_index, depth_index] +
    caa_table_12_1_df[i, "total_pax_this_period"]
}
dimnames( airport_countries_pax_array ) <- list( uk_airports, foreign_countries, seq(1,12,1))

# airport_regions_array
for (i in 1:nrow( caa_table_12_1_df ) ) {
  row_index <- match( caa_table_12_1_df[i, "UK_airport"], uk_airports )
  col_index <- match( caa_table_12_1_df[i, "foreign_region"], foreign_regions )
  depth_index <- as.numeric( stringr::str_sub( ( caa_table_12_1_df[ i, "this_period"] ), -2) )
  airport_regions_pax_array[row_index, col_index, depth_index] <- airport_regions_pax_array[row_index, col_index, depth_index] +
    caa_table_12_1_df[i, "total_pax_this_period"]
}
dimnames( airport_regions_pax_array ) <- list( uk_airports, foreign_regions, seq(1,12,1))

# Divide all passenger numbers by 2 as they include arrivals and departures
# Results in some non-integer values so clearly not as simple as a doubling of the count
airport_regions_pax_array <- airport_regions_pax_array / 2
airport_countries_pax_array <- airport_countries_pax_array / 2
airport_airport_pax_array <- airport_airport_pax_array / 2

# Check
HEATHROW <- sum( airport_airport_pax_array[ "HEATHROW",,])     #[1] 39414917
MANCHESTER <- sum( airport_airport_pax_array[ "MANCHESTER",,]) #[1] 14427330
HEATHROW/MANCHESTER #[1] 2.731962

saveRDS( airport_regions_pax_array
         , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_2024_foreign_region_to_uk_airport.rds" )
saveRDS( airport_airport_pax_array
         , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_2024_foreign_airport_to_uk_airport.rds" )
saveRDS( airport_countries_pax_array
         , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_2024_foreign_countries_to_uk_airport.rds" )

# Convert to % of total passenger numbers by month
percent_convert <- function( pax_array ){
  for(i in 1:12){
    pax_array[,,i] <- pax_array[,,i] / sum(pax_array[,,i] )
  }
}

# Make lookup table for UK airports and UTLA, LTLA, UK Region
#write.csv( uk_airports , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_2024_uk_airport.csv" )
# airport postcodes added manually
UK_airport_postcodes <- read.csv( "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_2024_uk_airport.csv" )
# ONS postcode to various geo regions/areas lookups (note these are only for England and Wales)
ONS_geo_lookups <- readr::read_csv( unzip( "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/pcd_oa_lsoa_msoa_ltla_utla_rgn_ctry_ew_may_2021_lu.zip"
                                    , "pcd_oa_lsoa_msoa_ltla_utla_rgn_ctry_ew_may_2021_lu_v2.csv"))
test <- subset( ONS_geo_lookups, ONS_geo_lookups$pcd %in% UK_airport_postcodes$postcode )
# Remove spaces from postcodes as creates difficulty in merging
remove_spaces <- function(x) gsub("\\s+", "", x)
UK_airport_postcodes$postcode <- remove_spaces( UK_airport_postcodes$postcode )
ONS_geo_lookups$pcd <- remove_spaces( ONS_geo_lookups$pcd )
Eng_Wales_airport_location <- merge( UK_airport_postcodes, ONS_geo_lookups
                              , by.x = "postcode", by.y = "pcd", all = FALSE)
saveRDS( Eng_Wales_airport_location , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/Eng_Wales_airport_location.rds" )


# 20 March 2025
# EV: For border crossings, can you work out how many arrivals we get per unit time in different ITL2 regions?
View( as.data.frame(airport_regions_pax_array[,,1] ))
View( as.data.frame(airport_regions_pax_array[,,2] ))

# Sum by airport for each month (array depth)
  # Check
  rowSums(airport_regions_pax_array[,,1])[1] == sum(airport_regions_pax_array[1,,1]) # [1] TRUE

  # Initialise new dataframe
  airport_pax_month <- data.frame( matrix( data = NA
                                          , nrow=nrow(airport_regions_pax_array)
                                          , ncol=12, 
                                          )
                                   )
  for(i in 1:12){
    airport_pax_month[ , i ] <- rowSums( airport_regions_pax_array[ , , i ] )
  }
  rownames( airport_pax_month ) <- rownames( airport_regions_pax_array )

# Check
  HEATHROW <- sum( airport_pax_month[ "HEATHROW",])     #[1] 39414917
  MANCHESTER <- sum( airport_pax_month[ "MANCHESTER",]) #[1] 14427330
  HEATHROW/MANCHESTER #[1] 2.731962
  
  
# then sum by month to get total arrivals by airport
  # Initialise new dataframe
  airport_pax_2024 <- as.data.frame( matrix( data = c( rownames( airport_pax_month ), rowSums( airport_pax_month ) )
                                             , nrow = nrow(airport_pax_month), ncol = 2
                                             )
  )
  colnames(airport_pax_2024) <- c("CAA_airport","pax_2024")
  
  # check
  sum( airport_regions_pax_array[1,,] ) == airport_pax_2024[1,2] # [1] TRUE
  
  HEATHROW <- subset( airport_pax_2024, CAA_airport == "HEATHROW" )$pax_2024    #[1] "39414916.5"
  MANCHESTER <- subset( airport_pax_2024, CAA_airport == "MANCHESTER" )$pax_2024 #[1] "14427329.5"
  as.numeric(HEATHROW) / as.numeric(MANCHESTER) #[1] 2.731962
  
  
# Add airport postcode
  airport_pax_2024_geo <- merge( x = Eng_Wales_airport_location, y = airport_pax_2024
                                , by.x = "CAA_airport", by.y = "CAA_airport"
                                , all.x = TRUE, all.y = TRUE
                                 )
  View( airport_pax_2024_geo )
  # Remove airports not in England & Wales (i.e. Scotland, Jersey, Guernsey etc)
  airport_pax_2024_geo_ew <- subset( airport_pax_2024_geo, !is.na( airport_pax_2024_geo$nat22cd ) )
  View( airport_pax_2024_geo_ew )

# Find ITL2
  LAD_to_ITL2 <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD_(December_2024)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2025)_Lookup_in_the_UK.csv")
  View(LAD_to_ITL2)
  airport_pax_2024_geo_ew_ITL2 <- merge(   x = airport_pax_2024_geo_ew
                                         , y = unique( LAD_to_ITL2[ , c("LAD24CD","LAD24NM","ITL225CD","ITL225NM")])
                                           , by.x = "ltla22cd", by.y = "LAD24CD"
                                           , all.x = TRUE , all.y = FALSE
                                           )
  View( airport_pax_2024_geo_ew_ITL2 )
# Sum passengers (pax) by ITL2
  length( unique( airport_pax_2024_geo_ew_ITL2$ITL225CD ) ) == nrow( airport_pax_2024_geo_ew_ITL2 ) # [1] FALSE
  table( airport_pax_2024_geo_ew_ITL2$ITL225CD )
  
  library(dplyr)
  # Collapse and sum using dplyr
  pax_2024_ew_ITL2 <- airport_pax_2024_geo_ew_ITL2 %>%
                        group_by( ITL225CD , ITL225NM ) %>%
                        summarise( pax_2024_total = sum( as.numeric( pax_2024 ) ) )
  pax_2024_ew_ITL2 <- pax_2024_ew_ITL2[ order( pax_2024_ew_ITL2[ , 1 ] ) , ]
  pax_2024_ew_ITL2[,"pax_2024_per_day"] <- pax_2024_ew_ITL2[,"pax_2024_total"] / 365
  pax_2024_ew_ITL2[,"pax_2024_per_day_integer"] <- round( pax_2024_ew_ITL2[,"pax_2024_per_day"] , digits = 0 )
  
  View(pax_2024_ew_ITL2)
  
  saveRDS( as.data.frame( pax_2024_ew_ITL2 ) , "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_pax_2024_ITL2.rds")
  test <- readRDS( "~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/CAA_pax_2024_ITL2.rds" ) 
  View(test)
  
  # Check
  sum( test$pax_2024_per_day_integer ) * 365 # [1] 117,744,620 # More because of rounding
  sum( test$pax_2024_total ) # [1] 117,744,262
  sum( airport_regions_pax_array[ rownames( airport_regions_pax_array ) %in% airport_pax_2024_geo_ew_ITL2$CAA_airport,,] ) # [1] 117,744,262
  
  # Check
  # Outer London - West and North West (TLI7) vs Greater Manchester (TLD3)
  OLWNW <- subset( pax_2024_ew_ITL2, pax_2024_ew_ITL2$ITL225CD == "TLI7")$pax_2024_per_day_integer
  #[1] 107986
  GM <- subset( pax_2024_ew_ITL2, pax_2024_ew_ITL2$ITL225CD == "TLD3")$pax_2024_per_day_integer
  #[1] 39527
  OLWNW / GM
  #[1] 2.731955
  HEATHROW    <- sum( subset( caa_table_12_1_df, UK_airport == "HEATHROW" )$total_pax_this_period ) / 2
  #[1] 39414917
  MANCHESTER <- sum( subset( caa_table_12_1_df, UK_airport == "MANCHESTER" )$total_pax_this_period ) / 2
  #[1] 14427330
  HEATHROW/MANCHESTER
  #[1] 2.731962
  

    # Plot map
  
  # Load required libraries
  library(ggplot2)
  library(sf)
  library(dplyr)
  
  # Load the shapefile for ITL2 regions (replace 'path_to_shapefile' with the actual file path)
  # https://ckan.publishing.service.gov.uk/dataset/international-territorial-level-2-january-2025-boundaries-uk-bgc/resource/78abd24e-5f15-4e04-a118-ddc649e9c769
  itl2_shapefile <- st_read("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2 shape files/ITL2_JAN_2025_UK_BFE.shp")
  
  # Example data frame with ITL2 region codes and values to color by
  region_data <- data.frame(
    ITL225CD = pax_2024_ew_ITL2$ITL225CD #c("UKC1", "UKC2", "UKD1", "UKD3"), # Replace with actual ITL2 codes
    , pax_2024_per_day_integer = pax_2024_ew_ITL2$pax_2024_per_day_integer #c(100, 200, 150, 300) # Replace with your column values
  )
  
  # Join the shapefile data with the region data
  itl2_map_data <- itl2_shapefile %>%
    left_join( region_data, by = c("ITL225CD" = "ITL225CD")) #c("ITL2_Code" = "ITL2_Code"))
  
  # Plot the map using ggplot2
  ggplot(data = itl2_map_data) +
    geom_sf(aes(fill = pax_2024_per_day_integer), color = "black") + # Fill regions by 'pax_2024_per_day_integer' column
    scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "white") + # Adjust color gradient
    theme_minimal() +
    labs(
      title = "UK Map with ITL2 Regions",
      fill = "International air passenger arrivals"
    )
  