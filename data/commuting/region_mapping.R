# Mapping from UTLA to ITL2

# ITL2
# https://geoportal.statistics.gov.uk/datasets/f2f482dc109e454599829bc5ab2c4418_0/explore
# https://geoportal.statistics.gov.uk/datasets/ons::international-territorial-level-2-january-2025-boundaries-uk-bfc/about


# Possible source for mapping?
# https://geoportal.statistics.gov.uk/datasets/15bdb6b0ff1c4a64b34a64e2e39f8caf_0/explore

# 2024 2025
ONS_geo_mapping_2025 <- read.csv( "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD_(December_2024)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2025)_Lookup_in_the_UK.csv")

# Check if LAU125CD = LAD24CD
length( ONS_geo_mapping_2025[,"LAU125CD"] ) # 362
sum( ONS_geo_mapping_2025[,"LAU125CD"] == ONS_geo_mapping_2025[,"LAD24CD"] ) # 329
sum( ONS_geo_mapping_2025[,"LAU125NM"] == ONS_geo_mapping_2025[,"LAD24NM"] ) # 360

# Differences are all in the codes for areas in Scotland
subset( ONS_geo_mapping_2025,
        ONS_geo_mapping_2025[,"LAU125CD"] != ONS_geo_mapping_2025[,"LAD24CD"])
subset( ONS_geo_mapping_2025,
        ONS_geo_mapping_2025[,"LAU125NM"] != ONS_geo_mapping_2025[,"LAD24NM"]) # Two differences, both in Scotland

# Conclusion: Ok to use either LAU125 or LAD24 if not using for areas in Scotland
# However, the ONS census 2021 commuting data has Scotland as a workplace so need to include

# Trim rows to England and Wales
#ITL2_LTLA_map_EW_2025 <- subset( ONS_geo_mapping_2025 , ONS_geo_mapping_2025$ITL125NM != "Scotland" )
# Trim columns to ITL2 and LAD
#ITL2_LTLA_map_EW_2025 <- ITL2_LTLA_map_EW_2025[,c("ITL225CD","ITL225NM","LAD24CD","LAD24NM")]
#saveRDS( ITL2_LTLA_map_EW_2025, "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2_LTLA_map_EW_2025.rds")

# Note that LAU121 (2021) and ITL321 (2021) area columns should be removed because
# these are smaller than the LTLAs for the LTLA = "Highland" and LAU125 = "Highland"
# and so there would be mutliple entries for the same LTLA leading to double counting of individuals
colnames( ONS_geo_mapping_2025 )
#[1] "ITL125CD" "ITL125NM" "ITL225CD" "ITL225NM" "ITL325CD" "ITL325NM" "LAU125CD" "LAU125NM" "LAD24CD"  "LAD24NM"
#[11] "ObjectId"
ONS_geo_mapping_2025_subset <- ONS_geo_mapping_2025[,c("LAD24CD","LAD24NM"
                                                       ,"ITL225CD","ITL225NM"
                                                       ,"ITL125CD","ITL125NM")]

# 2020 2021 may be more appropriate for ONS census data 2021
ONS_geo_mapping_2021 <- readxl::read_xlsx("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/LAD20_LAU121_ITL321_ITL221_ITL121_UK_LU_v2.xlsx")
ONS_geo_mapping_2021 <- as.data.frame( ONS_geo_mapping_2021 )

# Check if LAU125CD = LAD24CD
length( ONS_geo_mapping_2021[,"LAD20CD"] ) # 388
sum( ONS_geo_mapping_2021[,"LAU121CD"] == ONS_geo_mapping_2021[,"LAD20CD"] ) # 347
sum( ONS_geo_mapping_2021[,"LAU121NM"] == ONS_geo_mapping_2021[,"LAD20NM"] ) # 375

# Differences are all in the codes for areas in Scotland
subset( ONS_geo_mapping_2021,
        ONS_geo_mapping_2021[,"LAU121CD"] != ONS_geo_mapping_2021[,"LAD20CD"])
subset( ONS_geo_mapping_2021,
        ONS_geo_mapping_2021[,"LAU121NM"] != ONS_geo_mapping_2021[,"LAD20NM"]) # All differences are in Scotland

# Conclusion: Ok to use either LAU121 or LAD20 if not using for areas in Scotland
# However, the ONS census 2021 commuting data has Scotland as a workplace so need to include

# Trim rows to England and Wales
#ITL2_LTLA_map_EW_2021 <- subset( ONS_geo_mapping_2021 , ONS_geo_mapping_2021$ITL121NM != "Scotland" )
# Trim columns to ITL2 and LAD
#ITL2_LTLA_map_EW_2021 <- ITL2_LTLA_map_EW_2021[,c("ITL221CD","ITL221NM","LAD20CD","LAD20NM")]
#saveRDS( ITL2_LTLA_map_EW_2021, "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2_LTLA_map_EW_2021.rds")


# List of LTLAs used in ONS 2021 census data for commuting
ODWP01EW_LTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_LTLA.csv")
ODWP01EW_LTLA_residence <- unique( ODWP01EW_LTLA[,c("Lower.tier.local.authorities.code"
                                                    , "Lower.tier.local.authorities.label")] )
ODWP01EW_LTLA_workplace <- unique( ODWP01EW_LTLA[ , c("LTLA.of.workplace.code"
                                                      , "LTLA.of.workplace.label") ] )

# Sort by codes
ODWP01EW_LTLA_residence <- ODWP01EW_LTLA_residence[ order( ODWP01EW_LTLA_residence$Lower.tier.local.authorities.code ) , ]
ODWP01EW_LTLA_workplace <- ODWP01EW_LTLA_workplace[ order( ODWP01EW_LTLA_workplace$LTLA.of.workplace.code ) , ]

# Note that LAU121 (2021) and ITL321 (2021) area columns should be removed because
# these are smaller than the LTLAs for the LTLA = "Highland" and LAU125 = "Highland"
# and so there would be mutliple entries for the same LTLA leading to double counting of individuals
colnames( ONS_geo_mapping_2021 )
#[1] "LAD20CD"  "LAD20NM"  "LAU121CD" "LAU121NM" "ITL321CD" "ITL321NM" "ITL221CD" "ITL221NM" "ITL121CD" "ITL121NM"
ONS_geo_mapping_2021_subset <- ONS_geo_mapping_2021[,c("LAD20CD","LAD20NM"
                                                       ,"ITL221CD","ITL221NM"
                                                       ,"ITL121CD","ITL121NM")]

# Add 2020/21 ITL2
# Residence
area_mapping_residence <- merge( x = ODWP01EW_LTLA_residence, y = ONS_geo_mapping_2021_subset
                                 , by.x = "Lower.tier.local.authorities.code"
                                 , by.y = "LAD20CD"
                                 , all.x = TRUE , all.y = FALSE
                                 )
# Ensure have the LAD20CD column as well as the LTLA code column
area_mapping_residence <- merge( x = area_mapping_residence, y = ONS_geo_mapping_2021_subset[,c("LAD20NM","LAD20CD")]
                                 , by.x = "LAD20NM"
                                 , by.y = "LAD20NM"
                                 , all.x = TRUE , all.y = FALSE
)
# Workplace
area_mapping_workplace <- merge( x = ODWP01EW_LTLA_workplace, y = ONS_geo_mapping_2021_subset
                                 , by.x = "LTLA.of.workplace.code"
                                 , by.y = "LAD20CD"
                                 , all.x = TRUE , all.y = FALSE
)
# Ensure have the LAD20CD column as well as the LTLA code column
area_mapping_workplace <- merge( x = area_mapping_workplace, y = ONS_geo_mapping_2021_subset[,c("LAD20NM","LAD20CD")]
                                 , by.x = "LAD20NM"
                                 , by.y = "LAD20NM"
                                 , all.x = TRUE , all.y = FALSE
)

# Add 2024/25 ITL2
# Residence
area_mapping_residence <- merge( x = area_mapping_residence, y = ONS_geo_mapping_2025_subset
                                 , by.x = "Lower.tier.local.authorities.code"
                                 , by.y = "LAD24CD"
                                 , all.x = TRUE , all.y = FALSE
)
# Ensure have the LAD20CD column as well as the LTLA code column
area_mapping_residence <- merge( x = area_mapping_residence, y = ONS_geo_mapping_2025_subset[,c("LAD24NM","LAD24CD")]
                                 , by.x = "LAD24NM"
                                 , by.y = "LAD24NM"
                                 , all.x = TRUE , all.y = FALSE
)
# Workplace
area_mapping_workplace <- merge( x = area_mapping_workplace, y = ONS_geo_mapping_2025_subset
                                 , by.x = "LTLA.of.workplace.code"
                                 , by.y = "LAD24CD"
                                 , all.x = TRUE , all.y = FALSE
)
# Ensure have the LAD20CD column as well as the LTLA code column
area_mapping_workplace <- merge( x = area_mapping_workplace, y = ONS_geo_mapping_2025_subset[,c("LAD24NM","LAD24CD")]
                                 , by.x = "LAD24NM"
                                 , by.y = "LAD24NM"
                                 , all.x = TRUE , all.y = FALSE
)


# Remove rows only containing NA
area_mapping_residence <- area_mapping_residence[ rowSums( is.na( area_mapping_residence ) ) < ncol( area_mapping_residence ) , ]
area_mapping_workplace <- area_mapping_workplace[ rowSums( is.na( area_mapping_workplace ) ) < ncol( area_mapping_workplace ) , ]


# Re-order columns
area_mapping_residence <- area_mapping_residence[, c( "Lower.tier.local.authorities.code", "Lower.tier.local.authorities.label"
                                                      , "LAD20CD"  , "LAD20NM"
                                                      , "LAD24CD"  , "LAD24NM"
                                                      #, "LAU121CD" , "LAU121NM"
                                                      #, "LAU125CD" , "LAU125NM"
                                                      #, "ITL321CD" , "ITL321NM"
                                                      #, "ITL325CD" , "ITL325NM"
                                                      , "ITL221CD" , "ITL221NM"
                                                      , "ITL225CD" , "ITL225NM"
                                                      , "ITL121CD" , "ITL121NM"
                                                      , "ITL125CD" , "ITL125NM"
)]
area_mapping_workplace <- area_mapping_workplace[, c( "LTLA.of.workplace.code", "LTLA.of.workplace.label"
                                                      , "LAD20CD"  , "LAD20NM"
                                                      , "LAD24CD"  , "LAD24NM"
                                                      #, "LAU121CD" , "LAU121NM"
                                                      #, "LAU125CD" , "LAU125NM"
                                                      #, "ITL321CD" , "ITL321NM"
                                                      #, "ITL325CD" , "ITL325NM"
                                                      , "ITL221CD" , "ITL221NM"
                                                      , "ITL225CD" , "ITL225NM"
                                                      , "ITL121CD" , "ITL121NM"
                                                      , "ITL125CD" , "ITL125NM"
)]

saveRDS( area_mapping_residence , "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/area_mapping_ONS_commuting_residence.csv")
saveRDS( area_mapping_workplace , "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/area_mapping_ONS_commuting_workplace.csv")

# Looking at ITL2 changes between 2021 and 2025

ITL2_2021_to_2025 <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL_Level_2_(2021)_to_ITL_Level_2_(2025)_Lookup_in_the_UK.csv")


# Looking at another mapping file
# from https://open-geography-portalx-ons.hub.arcgis.com/datasets/ons::postcode-may-2021-to-oa-2021-to-lsoa-to-msoa-to-ltla-to-utla-to-rgn-to-ctry-may-2021-best-
ONS_mapping_info <- read.csv("~/mSCAPE/1_epidemic_modelling/Air traffic/Mar_2025_importation_prob/pcd_oa_lsoa_msoa_ltla_utla_rgn_ctry_ew_may_2021_lu_v2.csv")
sort( unique( ONS_mapping_info$msoa21cd ) )
