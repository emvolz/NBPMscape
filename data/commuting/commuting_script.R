# Inspecting the origin-destination data from the UK 2021 census,
# published by the Office for National Statistics at https://www.nomisweb.co.uk/sources/census_2021_od

# Origin-destination Workplace data
# These datasets show commuting flows between usual residence and place of work
# for people aged 16 years and over in employment or temporarily away from work
# in the week before Census Day. These data include one main flow and seven
# univariate datasets.

# LTLA
ODWP01EW_LTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_LTLA.csv")
length(table(ODWP01EW_LTLA$Lower.tier.local.authorities.label)) # 331 LTLAs

# UTLA
ODWP01EW_UTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_UTLA.csv")
length(table(ODWP01EW_UTLA$Upper.tier.local.authorities.label )) # 174 UTLAs

# Region
ODWP01EW_RGN <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_RGN.csv")
length(table(ODWP01EW_RGN$Regions.label)) # 10 regions
sum(ODWP01EW_RGN$Count) # 48,566,381 people
table(ODWP01EW_RGN$Place.of.work.indicator..4.categories..label)

# Create migration matrix for residence (rows) to workplace (columns)
# by Region
row_labels <- unique( ODWP01EW_RGN$Regions.label )
col_labels <- unique( ODWP01EW_RGN$Region.of.workplace.label )

commuting_Region_matrix <- matrix(0, nrow = length(row_labels), ncol = length(col_labels),
                                  dimnames = list(row_labels, col_labels))

for (i in seq_len( nrow( ODWP01EW_RGN ) ) ) {
  row <- ODWP01EW_RGN$Regions.label[ i ]
  col <- ODWP01EW_RGN$Region.of.workplace.label[ i ]
  value <- ODWP01EW_RGN$Count[i]

  commuting_Region_matrix[row, col] <- commuting_Region_matrix[row, col] + value
}

# Checks
sum(result_matrix[,'North East']) # 1,111,200
sum( subset(ODWP01EW_RGN, ODWP01EW_RGN$Region.of.workplace.label == "North East")$Count ) # 1,111,200
sum(result_matrix) # 48,566,381
sum(result_matrix[,'Does not apply']) # 20,792,708

# Invert matrix so matches the phylo prob matrix
commuting_Region_matrix <- t(commuting_Region_matrix)

saveRDS( commuting_Region_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_Region_matrix.rds" )

# by UTLA
row_labels <- unique( ODWP01EW_UTLA$Upper.tier.local.authorities.label )
col_labels <- unique( ODWP01EW_UTLA$UTLA.of.workplace.label )

commuting_UTLA_matrix <- matrix(0, nrow = length(row_labels), ncol = length(col_labels),
                               dimnames = list(row_labels, col_labels))

for (i in seq_len( nrow( ODWP01EW_UTLA ) ) ) {
  row <- ODWP01EW_UTLA$Upper.tier.local.authorities.label[ i ]
  col <- ODWP01EW_UTLA$UTLA.of.workplace.label[ i ]
  value <- ODWP01EW_UTLA$Count[i]

  commuting_UTLA_matrix[row, col] <- commuting_UTLA_matrix[row, col] + value
}

# Checks
sum(commuting_UTLA_matrix[,'Hartlepool']) # 33121
sum( subset(ODWP01EW_UTLA, ODWP01EW_UTLA$UTLA.of.workplace.label == "Hartlepool")$Count ) # 33121
sum(commuting_UTLA_matrix) # 48,566386
sum(commuting_UTLA_matrix[,'Does not apply']) # 20,792,708

# Invert matrix so matches the phylo prob matrix
commuting_UTLA_matrix <- t(commuting_UTLA_matrix)

saveRDS( commuting_UTLA_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_UTLA_matrix.rds" )

# by LTLA
row_labels <- unique( ODWP01EW_LTLA$Lower.tier.local.authorities.label )
col_labels <- unique( ODWP01EW_LTLA$LTLA.of.workplace.label )

commuting_LTLA_matrix <- matrix(0, nrow = length(row_labels), ncol = length(col_labels),
                                dimnames = list(row_labels, col_labels))

for (i in seq_len( nrow( ODWP01EW_LTLA ) ) ) {
  row <- ODWP01EW_LTLA$Lower.tier.local.authorities.label[ i ]
  col <- ODWP01EW_LTLA$LTLA.of.workplace.label[ i ]
  value <- ODWP01EW_LTLA$Count[i]

  commuting_LTLA_matrix[row, col] <- commuting_LTLA_matrix[row, col] + value
}

# Checks
sum(commuting_LTLA_matrix[,'Ipswich']) # 66386
sum( subset(ODWP01EW_LTLA, ODWP01EW_LTLA$LTLA.of.workplace.label == "Ipswich")$Count ) # 66386
sum(commuting_LTLA_matrix) # 48,566232
sum(commuting_LTLA_matrix[,'Does not apply']) # 20,792,704

# Invert matrix so matches the phylo prob matrix
commuting_LTLA_matrix <- t(commuting_LTLA_matrix)

saveRDS( commuting_LTLA_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_LTLA_matrix.rds" )

# Convert to commuting probability
commuting_Region_prob_matrix <- commuting_Region_matrix / rowSums( commuting_Region_matrix )
commuting_LTLA_prob_matrix <- commuting_LTLA_matrix / rowSums( commuting_LTLA_matrix )
commuting_UTLA_prob_matrix <- commuting_UTLA_matrix / rowSums( commuting_UTLA_matrix )

#Checks
sum( commuting_Region_prob_matrix[1,] )
sum( commuting_UTLA_prob_matrix[10,] )
sum( commuting_LTLA_prob_matrix[36,] )

# Save prob matrices
saveRDS( commuting_Region_prob_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_Region_prob_matrix.rds" )
saveRDS( commuting_UTLA_prob_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_UTLA_prob_matrix.rds" )
saveRDS( commuting_LTLA_prob_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_LTLA_prob_matrix.rds" )

#####################################
# Disaggregated by age
# Only MSOA available
ODWP04EW_MSOA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp04ew/ODWP04EW_MSOA.csv")
table(ODWP04EW_MSOA$Middle.layer.Super.Output.Areas.label) # 7264 areas
table(ODWP04EW_MSOA$Middle.layer.Super.Output.Areas.code)
colnames(ODWP04EW_MSOA)
unique(ODWP04EW_MSOA$Age..E...8.categories..label)

######################################

# Read in mapping dataframe to convert from LTLA to ITL2
# see 'region_mapping.R'
# Note that the year of the data matters as there are changes to the codes and names between years
# Here we are using the 2020/21 data which hopefully matches the 2021 ONS census data area codes
ITL2_LTLA_2021_map_EW <- readRDS( "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2_LTLA_2021_map_EW.rds")
#ITL2_LTLA_map_EW <- readRDS("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2_LTLA_map_EW.rds")
#length( unique( ITL2_LTLA_map_EW$LAD24NM ) ) # 329 LTLAs
length( unique( ITL2_LTLA_2021_map_EW$LAD20NM ) ) # 347 LTLAs

# Read in LTLA matrices
commuting_LTLA_prob_matrix <- readRDS( "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_LTLA_prob_matrix.rds" )
commuting_LTLA_matrix      <- readRDS( "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_LTLA_matrix.rds" )

#sort( colnames(commuting_LTLA_matrix ) ) == unique( ITL2_LTLA_map$LAD24NM )
sort( colnames(commuting_LTLA_matrix ) ) == sort( unique( ITL2_LTLA_2021_map_EW$LAD20NM ) ) # not the same length

# Load ONS commuting data for LTLAs
ODWP01EW_LTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_LTLA.csv")
length( table( ODWP01EW_LTLA$Lower.tier.local.authorities.label ) ) # 331 LTLAs

# Merge to add ITL2 areas to commuting data for RESIDENCE
ODWP01EW_ITL2 <- merge( x = ODWP01EW_LTLA , y = ITL2_LTLA_2021_map_EW
                        , by.x = "Lower.tier.local.authorities.code", by.y = "LAD20CD"
                        , all.x = TRUE , all.y = FALSE)

# Checks
sum( ODWP01EW_ITL2$Count ) # 48566232
sum( ODWP01EW_LTLA$Count ) # 48566232
unique( ODWP01EW_LTLA$Lower.tier.local.authorities.code ) == unique( ODWP01EW_ITL2$Lower.tier.local.authorities.code )
sort( unique( ODWP01EW_ITL2$Lower.tier.local.authorities.label ) ) == sort( unique( ODWP01EW_ITL2$LAD20NM ) ) # Not all match
subset( ODWP01EW_ITL2, is.na( ODWP01EW_ITL2$LAD20NM ) ) #

# If checks ok then rename columns
names(ODWP01EW_ITL2)[8] <- "ITL221CD_residence"
names(ODWP01EW_ITL2)[9] <- "ITL221NM_residence"
names(ODWP01EW_ITL2)[10] <- "LAD20NM_residence"

# Merge to add ITL2 areas to commuting data for WORKPLACE
ODWP01EW_ITL2 <- merge( x = ODWP01EW_ITL2 , y = ITL2_LTLA_2021_map_EW
                        , by.x = "LTLA.of.workplace.code", by.y = "LAD20CD"
                        , all.x = TRUE , all.y = FALSE)

# Checks
sum( ODWP01EW_ITL2$Count ) # 48566232
sum( ODWP01EW_LTLA$Count ) # 48566232
sort( unique( ODWP01EW_LTLA$Lower.tier.local.authorities.code ) ) == sort( unique( ODWP01EW_ITL2$Lower.tier.local.authorities.code ) ) # All true
sort( unique( ODWP01EW_ITL2$Lower.tier.local.authorities.label ) ) == sort( unique( ODWP01EW_ITL2$LAD20NM ) ) # Not the same length
subset( ODWP01EW_ITL2, is.na( ODWP01EW_ITL2$LAD20NM ) ) # Lots of NA rows

# Rename columns
names(ODWP01EW_ITL2)[11] <- "ITL221CD_workplace"
names(ODWP01EW_ITL2)[12] <- "ITL221NM_workplace"
names(ODWP01EW_ITL2)[13] <- "LAD20NM_workplace"

# Adjust ITL2 workplace for LTLA = "Does not apply"
for ( i in 1:nrow( ODWP01EW_ITL2 ) ){
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == -8 &
      is.na( ODWP01EW_ITL2$ITL221CD_workplace[i] ) &
      is.na( ODWP01EW_ITL2$ITL221NM_workplace[i] ) &
      is.na( ODWP01EW_ITL2$LAD20NM_workplace[i] ) ){

    ODWP01EW_ITL2$ITL221CD_workplace[i] = -8
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Does not apply"
    ODWP01EW_ITL2$LAD20NM_workplace[i] = "Does not apply"
  }
}

# Adjust ITL2 workplace for LTLA = "Workplace is offshore installation"
for ( i in 1:nrow( ODWP01EW_ITL2 ) ){
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == 888888888 &
      is.na( ODWP01EW_ITL2$ITL221CD_workplace[i] ) &
      is.na( ODWP01EW_ITL2$ITL221NM_workplace[i] ) &
      is.na( ODWP01EW_ITL2$LAD20NM_workplace[i] ) ){

    ODWP01EW_ITL2$ITL221CD_workplace[i] = 888888888
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Workplace is offshore installation"
    ODWP01EW_ITL2$LAD20NM_workplace[i] = "Workplace is offshore installation"
  }
}

# Adjust ITL2 workplace for LTLA = "Workplace is outside the UK"
for ( i in 1:nrow( ODWP01EW_ITL2 ) ){
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == 999999999 &
      is.na( ODWP01EW_ITL2$ITL221CD_workplace[i] ) &
      is.na( ODWP01EW_ITL2$ITL221NM_workplace[i] ) &
      is.na( ODWP01EW_ITL2$LAD20NM_workplace[i] ) ){

    ODWP01EW_ITL2$ITL221CD_workplace[i] = 999999999
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Workplace is outside the UK"
    ODWP01EW_ITL2$LAD20NM_workplace[i] = "Workplace is outside the UK"
  }
}

# Checking LTLAs for which there is no ITL2 conversion
missing_ITL2 <- subset( ODWP01EW_ITL2, is.na( ODWP01EW_ITL2$ITL221CD_residence ) )
LTLA_with_no_ITL2 <- as.data.frame( unique( missing_ITL2$Lower.tier.local.authorities.code ) )
names(LTLA_with_no_ITL2) = "Lower.tier.local.authorities.code"
LTLA_with_no_ITL2 <- merge( x = LTLA_with_no_ITL2 , y = missing_ITL2[,c("Lower.tier.local.authorities.code","Lower.tier.local.authorities.label")]
                            , by.x = "Lower.tier.local.authorities.code", by.y = "Lower.tier.local.authorities.code"
                            , all.x = FALSE, all.y = FALSE )
LTLA_with_no_ITL2 <- unique( LTLA_with_no_ITL2[, c(1,2)])

extras <- merge( x = LTLA_with_no_ITL2, y = ITL2_LTLA_map
                , by.x = "Lower.tier.local.authorities.code" , by.y = "LAD24CD"
                , all.x = TRUE, all.y = FALSE  )
# E06000061 North Northamptonshire and E06000062 West Northamptonshire changed between 2020 and 2024


# Merge in missing geo areas
for( i in 1:nrow( ODWP01EW_ITL2 ) ){

  # residence
  if( ODWP01EW_ITL2$Lower.tier.local.authorities.code[i] %in% extras$Lower.tier.local.authorities.code ){
    extras_subset_residence = subset( extras
                            , extras$Lower.tier.local.authorities.code == ODWP01EW_ITL2$Lower.tier.local.authorities.code[i] )
    ODWP01EW_ITL2$ITL221CD_residence[i] = extras_subset_residence$ITL225CD
    ODWP01EW_ITL2$ITL221NM_residence[i] = extras_subset_residence$ITL225NM
  }

  # workplace
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] %in% extras$Lower.tier.local.authorities.code ){
    extras_subset_workplace = subset( extras
                            , extras$Lower.tier.local.authorities.code == ODWP01EW_ITL2$LTLA.of.workplace.code[i] )
    ODWP01EW_ITL2$ITL221CD_workplace[i] = extras_subset_workplace$ITL225CD
    ODWP01EW_ITL2$ITL221NM_workplace[i] = extras_subset_workplace$ITL225NM
  }

}

# Checking LTLAs for which there is no ITL2 conversion
missing_ITL2_after <- subset( ODWP01EW_ITL2, is.na( ODWP01EW_ITL2$ITL221CD_residence ) )

# Make commuting matrix for ITL2
row_labels <- unique( ODWP01EW_ITL2$ITL221CD_residence )
col_labels <- unique( ODWP01EW_ITL2$ITL221CD_workplace )

commuting_ITL2_2021cd_matrix <- matrix(0, nrow = length(row_labels), ncol = length(col_labels),
                                       dimnames = list(row_labels, col_labels))

for (i in seq_len( nrow( ODWP01EW_ITL2 ) ) ) {
  print(i)
  row <- ODWP01EW_ITL2$ITL221CD_residence[ i ]
  col <- ODWP01EW_ITL2$ITL221CD_workplace[ i ]
  value <- ODWP01EW_ITL2$Count[i]

  commuting_ITL2_2021cd_matrix[row, col] <- commuting_ITL2_2021cd_matrix[row, col] + value
}

# Checks
sum(commuting_LTLA_matrix[,'Ipswich']) # 66386
sum( subset(ODWP01EW_LTLA, ODWP01EW_LTLA$LTLA.of.workplace.label == "Ipswich")$Count ) # 66386
sum(commuting_LTLA_matrix) # 48,566232
sum(commuting_LTLA_matrix[,'Does not apply']) # 20,792,704

# Invert matrix so matches the phylo prob matrix
commuting_LTLA_matrix <- t(commuting_LTLA_matrix)

saveRDS( commuting_LTLA_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_LTLA_matrix.rds" )

# Convert to commuting probability
commuting_Region_prob_matrix <- commuting_Region_matrix / rowSums( commuting_Region_matrix )

# Data now has ITL2 from 2021 but now also add ITL2 classification for 2025
ITL2_2021_to_2025 <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL_Level_2_(2021)_to_ITL_Level_2_(2025)_Lookup_in_the_UK.csv")

################################################################################

# Read in mapping for ONS commuting data residence and workplace LTLA to other geo area categories (and for different years)
# see 'region_mapping.R'
area_mapping_residence <- readRDS( "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/area_mapping_ONS_commuting_residence.csv")
area_mapping_workplace <- readRDS( "~/mSCAPE/1_epidemic_modelling/UK_region_lookups/area_mapping_ONS_commuting_workplace.csv")

length(unique(area_mapping_residence$Lower.tier.local.authorities.code)) # 331
length(unique(area_mapping_workplace$LTLA.of.workplace.code)) # 377

nrow(area_mapping_residence) == length(unique(area_mapping_residence$Lower.tier.local.authorities.code)) # TRUE
nrow(area_mapping_workplace) == length(unique(area_mapping_workplace$LTLA.of.workplace.code)) # FALSE

# Note that LAU121 (2021) and ITL321 (2021) area columns should be removed because
# these are smaller than the LTLAs for the LTLA = "Highland" and LAU125 = "Highland"
# and so there would be mutliple entries for the same LTLA leading to double counting of individuals
area_mapping_residence_subset <- area_mapping_residence[,c("Lower.tier.local.authorities.code","Lower.tier.local.authorities.label"
                                                           ,"LAD20CD","LAD20NM","LAD24CD","LAD24NM"
                                                           #,"LAU121CD","LAU121NM"
                                                           #,"LAU125CD","LAU125NM"
                                                           ,"ITL221CD","ITL221NM","ITL225CD","ITL225NM")]
area_mapping_workplace_subset <- area_mapping_workplace[,c("LTLA.of.workplace.code","LTLA.of.workplace.label"
                                                           ,"LAD20CD","LAD20NM","LAD24CD","LAD24NM"
                                                           #,"LAU121CD","LAU121NM"
                                                           #,"LAU125CD","LAU125NM"
                                                           ,"ITL221CD","ITL221NM","ITL225CD","ITL225NM")]

colnames( area_mapping_residence_subset ) <- paste( colnames( area_mapping_residence_subset )
                                                   , 'residence', sep = '_' )
colnames( area_mapping_workplace_subset ) <- paste( colnames( area_mapping_workplace_subset )
                                                    , 'workplace', sep = '_' )


nrow(area_mapping_residence_subset) == length(unique(area_mapping_residence$Lower.tier.local.authorities.code)) # TRUE
nrow(area_mapping_workplace_subset) == length(unique(area_mapping_workplace$LTLA.of.workplace.code)) # FALSE

View(table(area_mapping_workplace_subset$LTLA.of.workplace.code_workplace))
t_LTLA <- table(area_mapping_workplace_subset$LTLA.of.workplace.code_workplace)

multiple_entry_LTLA <- t_LTLA[ t_LTLA > 1]

# Load ONS commuting data for LTLAs
ODWP01EW_LTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_LTLA.csv")

ODWP01EW_ITL2 <- merge( x = ODWP01EW_LTLA , y = area_mapping_residence_subset
                          , by.x = "Lower.tier.local.authorities.code", by.y ="Lower.tier.local.authorities.code_residence"
                          , all.x = TRUE, all.y = FALSE
                        )

ODWP01EW_ITL2 <- merge( x = ODWP01EW_ITL2 , y = area_mapping_workplace_subset
                        , by.x = "LTLA.of.workplace.code", by.y ="LTLA.of.workplace.code_workplace"
                        , all.x = TRUE, all.y = FALSE
)

# Adjust for rows where the workplace is a non-geographic option
for ( i in 1:nrow( ODWP01EW_ITL2 ) ){
  # Adjust ITL2 workplace for LTLA = "Does not apply"
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == -8  ){
    ODWP01EW_ITL2$LAD20CD_workplace[i]  = -8
    ODWP01EW_ITL2$LAD24CD_workplace[i]  = -8
    ODWP01EW_ITL2$LAU121CD_workplace[i] = -8
    ODWP01EW_ITL2$LAU125CD_workplace[i] = -8
    ODWP01EW_ITL2$ITL221CD_workplace[i] = -8
    ODWP01EW_ITL2$ITL225CD_workplace[i] = -8

    ODWP01EW_ITL2$LAD20NM_workplace[i]  = "Does not apply"
    ODWP01EW_ITL2$LAD24NM_workplace[i]  = "Does not apply"
    ODWP01EW_ITL2$LAU121NM_workplace[i] = "Does not apply"
    ODWP01EW_ITL2$LAU125NM_workplace[i] = "Does not apply"
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Does not apply"
    ODWP01EW_ITL2$ITL225NM_workplace[i] = "Does not apply"
  }

  # Adjust ITL2 workplace for LTLA = "Workplace is offshore installation"
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == 888888888  ){
    ODWP01EW_ITL2$LAD20CD_workplace[i]  = 888888888
    ODWP01EW_ITL2$LAD24CD_workplace[i]  = 888888888
    ODWP01EW_ITL2$LAU121CD_workplace[i] = 888888888
    ODWP01EW_ITL2$LAU125CD_workplace[i] = 888888888
    ODWP01EW_ITL2$ITL221CD_workplace[i] = 888888888
    ODWP01EW_ITL2$ITL225CD_workplace[i] = 888888888

    ODWP01EW_ITL2$LAD20NM_workplace[i]  = "Workplace is offshore installation"
    ODWP01EW_ITL2$LAD24NM_workplace[i]  = "Workplace is offshore installation"
    ODWP01EW_ITL2$LAU121NM_workplace[i] = "Workplace is offshore installation"
    ODWP01EW_ITL2$LAU125NM_workplace[i] = "Workplace is offshore installation"
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Workplace is offshore installation"
    ODWP01EW_ITL2$ITL225NM_workplace[i] = "Workplace is offshore installation"
  }

  # Adjust ITL2 workplace for LTLA = "Workplace is outside the UK"
  if( ODWP01EW_ITL2$LTLA.of.workplace.code[i] == 999999999  ){
    ODWP01EW_ITL2$LAD20CD_workplace[i]  = 999999999
    ODWP01EW_ITL2$LAD24CD_workplace[i]  = 999999999
    ODWP01EW_ITL2$LAU121CD_workplace[i] = 999999999
    ODWP01EW_ITL2$LAU125CD_workplace[i] = 999999999
    ODWP01EW_ITL2$ITL221CD_workplace[i] = 999999999
    ODWP01EW_ITL2$ITL225CD_workplace[i] = 999999999

    ODWP01EW_ITL2$LAD20NM_workplace[i]  = "Workplace is outside the UK"
    ODWP01EW_ITL2$LAD24NM_workplace[i]  = "Workplace is outside the UK"
    ODWP01EW_ITL2$LAU121NM_workplace[i] = "Workplace is outside the UK"
    ODWP01EW_ITL2$LAU125NM_workplace[i] = "Workplace is outside the UK"
    ODWP01EW_ITL2$ITL221NM_workplace[i] = "Workplace is outside the UK"
    ODWP01EW_ITL2$ITL225NM_workplace[i] = "Workplace is outside the UK"
  }
}

# Inspect any rows containing NA in any column
# Remove rows only containing NA
View( ODWP01EW_ITL2[ rowSums( is.na( ODWP01EW_ITL2 ) ) > 0 , ] )

sum(ODWP01EW_ITL2$Count) # [1] 48,588,068 #** More than 21,000 counts added **
sum(ODWP01EW_LTLA$Count) # [1] 48,566,232

################################################################################
# Trying another way

# Load ONS commuting data for LTLAs
ODWP01EW_LTLA <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/2021/odwp/odwp/odwp01ew/ODWP01EW_LTLA.csv")

# from 'region_mapping.R'
ONS_geo_mapping_2021_subset
ONS_geo_mapping_2025_subset

# Add ITL2 columns (2021 and 2025) for residence and workplace
new_col_start <- ncol(ODWP01EW_LTLA) + 1
ODWP01EW_ITL2 <- ODWP01EW_LTLA
ODWP01EW_ITL2[, seq( new_col_start, new_col_start + 7, 1 ) ] <- NA
colnames(ODWP01EW_ITL2)[seq( new_col_start, new_col_start + 3, 1 )] <- c("ITL221CD_residence","ITL221NM_residence","ITL225CD_residence","ITL225NM_residence")
colnames(ODWP01EW_ITL2)[seq( new_col_start+4, new_col_start + 7, 1 )] <- c("ITL221CD_workplace","ITL221NM_workplace","ITL225CD_workplace","ITL225NM_workplace")

# Cycle through rows of dataframe
for (i in 1: nrow( ODWP01EW_ITL2 ) ){
  #print(i)
  # Residence info
  LTLA_CD_residence <- ODWP01EW_ITL2$Lower.tier.local.authorities.code[ i ]
  residence_matching_row_2021 <- subset( ONS_geo_mapping_2021_subset, ONS_geo_mapping_2021_subset$LAD20CD == LTLA_CD_residence )
  residence_matching_row_2025 <- subset( ONS_geo_mapping_2025_subset, ONS_geo_mapping_2025_subset$LAD24CD == LTLA_CD_residence )

  # Add ITL2 codes
  # If there is an ITL2 code then add it, otherwise use NA
  ODWP01EW_ITL2$ITL221CD_residence[ i ] <- ifelse( length( residence_matching_row_2021$ITL221CD ) == 0, NA, residence_matching_row_2021$ITL221CD )
  ODWP01EW_ITL2$ITL225CD_residence[ i ] <- ifelse( length( residence_matching_row_2025$ITL225CD ) == 0, NA, residence_matching_row_2025$ITL225CD )
  # Add ITL2 names
  ODWP01EW_ITL2$ITL221NM_residence[ i ] <- ifelse( length( residence_matching_row_2021$ITL221NM ) == 0, NA, residence_matching_row_2021$ITL221NM )
  ODWP01EW_ITL2$ITL225NM_residence[ i ] <- ifelse( length(residence_matching_row_2025$ITL225NM) == 0, NA, residence_matching_row_2025$ITL225NM )

  # Workplace codes
  LTLA_CD_workplace <- ODWP01EW_ITL2$LTLA.of.workplace.code[ i ]
  # Adjust ITL2 workplace for LTLA = "Does not apply"
  if( LTLA_CD_workplace == -8  ){
    ODWP01EW_ITL2[ i, c("ITL221CD_workplace","ITL225CD_workplace")]  = -8
    ODWP01EW_ITL2[ i, c("ITL221NM_workplace","ITL225NM_workplace")]  = "Does not apply"
  } else if ( LTLA_CD_workplace == 888888888 ){
    ODWP01EW_ITL2[ i, c("ITL221CD_workplace","ITL225CD_workplace")]  = 888888888
    ODWP01EW_ITL2[ i, c("ITL221NM_workplace","ITL225NM_workplace")]  = "Workplace is offshore installation"
  } else if ( LTLA_CD_workplace == 999999999 ){
    ODWP01EW_ITL2[ i, c("ITL221CD_workplace","ITL225CD_workplace")]  = 999999999
    ODWP01EW_ITL2[ i, c("ITL221NM_workplace","ITL225NM_workplace")]  = "Workplace is outside the UK"
  } else {
    workplace_matching_row_2021 <- subset( ONS_geo_mapping_2021_subset, ONS_geo_mapping_2021_subset$LAD20CD == LTLA_CD_workplace )
    workplace_matching_row_2025 <- subset( ONS_geo_mapping_2025_subset, ONS_geo_mapping_2025_subset$LAD24CD == LTLA_CD_workplace )
    # Add ITL2 codes
    ODWP01EW_ITL2$ITL221CD_workplace[ i ] <- ifelse( length( workplace_matching_row_2021$ITL221CD ) == 0, NA, workplace_matching_row_2021$ITL221CD )
    ODWP01EW_ITL2$ITL225CD_workplace[ i ] <- ifelse( length( workplace_matching_row_2025$ITL225CD ) == 0, NA, workplace_matching_row_2025$ITL225CD )
    # Add ITL2 names
    ODWP01EW_ITL2$ITL221NM_workplace[ i ] <- ifelse( length( workplace_matching_row_2021$ITL221NM ) == 0, NA, workplace_matching_row_2021$ITL221NM )
    ODWP01EW_ITL2$ITL225NM_workplace[ i ] <- ifelse( length( workplace_matching_row_2025$ITL225NM ) == 0, NA, workplace_matching_row_2025$ITL225NM )
  }
}

# Check
unique( ODWP01EW_ITL2$ITL221CD_residence )
unique( ODWP01EW_ITL2$ITL225CD_residence )

# View rows with ITL2 = NA
View( ODWP01EW_ITL2[ rowSums( is.na( ODWP01EW_ITL2 ) ) > 1 , ] )

# Some rows only have a 2021 ITL2 and some only have a 2025 ITL2
# So add a column to take one or the other
new_col_start <- ncol(ODWP01EW_LTLA) + 1
ODWP01EW_ITL2[, seq( new_col_start, new_col_start+3, 1 ) ] <- NA
colnames(ODWP01EW_ITL2)[seq( new_col_start, new_col_start+3, 1 )] <- c("ITL2_CD_residence_final","ITL2_NM_residence_final","ITL2_CD_workplace_final","ITL2_NM_workplace_final")

for (i in 1: nrow( ODWP01EW_ITL2 ) ){
  ODWP01EW_ITL2$ITL2_CD_residence_final[ i ] <- ifelse( is.na( ODWP01EW_ITL2$ITL225CD_residence[ i ] ), ODWP01EW_ITL2$ITL221CD_residence[ i ] , ODWP01EW_ITL2$ITL225CD_residence[ i ] )
  ODWP01EW_ITL2$ITL2_NM_residence_final[ i ] <- ifelse( is.na( ODWP01EW_ITL2$ITL225NM_residence[ i ] ), ODWP01EW_ITL2$ITL221NM_residence[ i ] , ODWP01EW_ITL2$ITL225NM_residence[ i ] )
  ODWP01EW_ITL2$ITL2_CD_workplace_final[ i ] <- ifelse( is.na( ODWP01EW_ITL2$ITL225CD_workplace[ i ] ), ODWP01EW_ITL2$ITL221CD_workplace[ i ] , ODWP01EW_ITL2$ITL225CD_workplace[ i ] )
  ODWP01EW_ITL2$ITL2_NM_workplace_final[ i ] <- ifelse( is.na( ODWP01EW_ITL2$ITL225NM_workplace[ i ] ), ODWP01EW_ITL2$ITL221NM_workplace[ i ] , ODWP01EW_ITL2$ITL225NM_workplace[ i ] )
}

# Checks
sort( unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) )
sort( unique( ODWP01EW_ITL2$ITL2_NM_residence_final ) )
sort( unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) )
sort( unique( ODWP01EW_ITL2$ITL2_NM_workplace_final ) )

# Some final ITL2 codes are not in the ITL2 2025 list
ITL2_2021_to_2025 <- read.csv("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL_Level_2_(2021)_to_ITL_Level_2_(2025)_Lookup_in_the_UK.csv")
unique( ITL2_2021_to_2025$ITL225CD )

# The changes from 2021 to 2025 in ITL2 codes are not one-to-one:
# - TLC1 is split between TLC3 and TLC4
# - TLH1 is split between TLH4, TLH5 and TLH6
# - TLK1 is split between TLK5, TLK6 and TLK7
# - TLL1 is split between TLL3, TLL4, and TLL5
# - TLL2 is split between TLL3, TLL4
# Note that the reverse direction is also not one-to-one

# Residence codes not in 2025 ITL2 list
'%ni%' <- Negate('%in%')
missing_ITL2_index <- which( sort( unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) ) %ni% sort( unique( ITL2_2021_to_2025$ITL225CD ) ) )
sort( unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) )[ missing_ITL2_index ] # TLK2

# Workplace codes not in 2025 ITL2 list
'%ni%' <- Negate('%in%')
missing_ITL2_index <- which( sort( unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) ) %ni% sort( unique( ITL2_2021_to_2025$ITL225CD ) ) )
sort( unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) )[ missing_ITL2_index ] #"-8"        "888888888" "999999999" "TLK2"

# From ITL2_2021_to_2025 we can see that TLK2 maps to TLK6
ITL2_update_row <- subset( ITL2_2021_to_2025, ITL2_2021_to_2025$ITL221CD == "TLK2" )

# Replace TLK2 with TLK6 in final ITL2 columns
for( i in 1:nrow( ODWP01EW_ITL2 ) ){
  # Residence
  # Code
  ODWP01EW_ITL2$ITL2_CD_residence_final[ i ] <- ifelse( ODWP01EW_ITL2$ITL2_CD_residence_final[ i ] == ITL2_update_row$ITL221CD , ITL2_update_row$ITL225CD,  ODWP01EW_ITL2$ITL2_CD_residence_final[ i ] )
  # Name
  ODWP01EW_ITL2$ITL2_NM_residence_final[ i ] <- ifelse( ODWP01EW_ITL2$ITL2_NM_residence_final[ i ] == ITL2_update_row$ITL221CD , ITL2_update_row$ITL225NM,  ODWP01EW_ITL2$ITL2_NM_residence_final[ i ] )
  # Workplace
  # Code
  ODWP01EW_ITL2$ITL2_CD_workplace_final[ i ] <- ifelse( ODWP01EW_ITL2$ITL2_CD_workplace_final[ i ] == ITL2_update_row$ITL221CD , ITL2_update_row$ITL225CD,  ODWP01EW_ITL2$ITL2_CD_workplace_final[ i ] )
  # Name
  ODWP01EW_ITL2$ITL2_NM_workplace_final[ i ] <- ifelse( ODWP01EW_ITL2$ITL2_NM_workplace_final[ i ] == ITL2_update_row$ITL221CD , ITL2_update_row$ITL225NM,  ODWP01EW_ITL2$ITL2_NM_workplace_final[ i ] )

}

# Checks
unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) %in% "TLK2"
unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) %in% "TLK2"

sort( unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) )
sort( unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) )

redundant_ITL2 <- ITL2_update_row$ITL221CD[ ITL2_update_row$ITL221CD %ni% ITL2_update_row$ITL225CD ]

# Finally, generate the commuting matrix using ITL2 geo areas
# by ITL2
row_labels <- sort( unique( ODWP01EW_ITL2$ITL2_CD_residence_final ) )
col_labels <- sort( unique( ODWP01EW_ITL2$ITL2_CD_workplace_final ) )

commuting_ITL2_matrix <- matrix(0, nrow = length(row_labels), ncol = length(col_labels),
                                 dimnames = list(row_labels, col_labels))

for (i in seq_len( nrow( ODWP01EW_ITL2 ) ) ) {
  row <- ODWP01EW_ITL2$ITL2_CD_residence_final[ i ]
  col <- ODWP01EW_ITL2$ITL2_CD_workplace_final[ i ]
  value <- ODWP01EW_ITL2$Count[i]

  commuting_ITL2_matrix[row, col] <- commuting_ITL2_matrix[row, col] + value
}

# Checks
sum(commuting_ITL2_matrix[,'TLK6']) # 713254
sum( subset( ODWP01EW_ITL2, ODWP01EW_ITL2$ITL2_CD_workplace_final == "TLK6")$Count ) # 713254
sum(commuting_ITL2_matrix) # 48,566232
sum(commuting_ITL2_matrix[,"-8"]) # 20,792,704

# Convert to commuting probability
commuting_ITL2_prob_matrix <- commuting_ITL2_matrix / rowSums( commuting_ITL2_matrix )
# Check
sum( commuting_ITL2_prob_matrix[1,]) # [1] 1 # Probabilities from the same origin have to sum to 1
sum( commuting_ITL2_prob_matrix[30,]) # [1] 1 # Probabilities from the same origin have to sum to 1

# Transpose matrix so matches the phylo prob matrix (rows = destination/workplace, cols = origin/residence)
commuting_ITL2_matrix      <- t(commuting_ITL2_matrix)
commuting_ITL2_prob_matrix <- t(commuting_ITL2_prob_matrix)

# Check
sum( commuting_ITL2_prob_matrix[,1]) # [1] 1 # Probabilities from the same origin have to sum to 1
sum( commuting_ITL2_prob_matrix[,30]) # [1] 1 # Probabilities from the same origin have to sum to 1

saveRDS( commuting_ITL2_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_ITL2_matrix.rds" )
saveRDS( commuting_ITL2_prob_matrix, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/commuting_ITL2_prob_matrix.rds" )

# Key for codes and names
ITL2_residence <- as.data.frame( unname( ODWP01EW_ITL2[ , c("ITL2_CD_residence_final", "ITL2_NM_residence_final" ) ] ) )
names(ITL2_residence) <- c("code","name")
ITL2_workplace <- as.data.frame( unname( ODWP01EW_ITL2[ , c("ITL2_CD_workplace_final", "ITL2_NM_workplace_final" ) ] ) )
names(ITL2_workplace) <- c("code","name")
ITL2_key <- data.frame(matrix(ncol = 2, nrow = 1))
names(ITL2_key) <- c("code","name")
ITL2_key <-  rbind( ITL2_key, ITL2_residence)
ITL2_key <-  rbind( ITL2_key, ITL2_workplace)
ITL2_key <- ITL2_key[-1,]
ITL2_key <-  unique( ITL2_key )
ITL2_key <-  ITL2_key[ order( ITL2_key$code ) , ]
row.names(ITL2_key) <- 1:nrow(ITL2_key)
saveRDS( ITL2_key, "~/mSCAPE/1_epidemic_modelling/UK_census_origin_destination/ITL2_key.rds" )


# Plot map
#** NOT COMPLETE**
# Load required libraries
library(ggplot2)
library(sf)
library(dplyr)

# Read in final commuting probability data
commuting_ITL2_inprob_list <- readRDS("~/GitHub/NBPMscape/data/commuting_ITL2_inprob_list.rds")
View(commuting_ITL2_inprob_list)

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


# Plot day of week contact rate probability
dowcont = c(0.1043502,0.1402675,0.1735913,0.1437642,0.1596205,0.1445298,0.1338766)
# Load ggplot2 library
library(ggplot2)

# Create sample data
dow_df <- data.frame(
  category = c("Sun", "Mon", "Tue", "Wed", "Thur", "Fri", "Sat"),
  value = dowcont
)

# Create the bar plot
ggplot(dow_df, aes(x = category, y = value)) +
  geom_bar(stat = "identity", fill = "blue") +
  scale_x_discrete(limits = dow_df$category)
  #theme_minimal() +
    theme(
      axis.text = element_text(size = 20),       # Increase axis text size
      axis.title = element_text(size = 20),      # Increase axis title size
      plot.title = element_text(size = 18)       # Increase plot title size
    ) +
  labs(title = "Basic Bar Plot"
       ,x = "Category"
       ,y = "Probaility"
       )



# Read in final commuting probability data
commuting_ITL2_inprob_list <- readRDS("~/GitHub/NBPMscape/data/commuting_ITL2_inprob_list.rds")
View(commuting_ITL2_inprob_list)
commuting_ITL2_inprob_matrix <- do.call(cbind, commuting_ITL2_inprob_list)
commuting_ITL2_inprob_df <- as.data.frame( commuting_ITL2_inprob_matrix )

library(sf)
library(ggplot2)
library(tidyr)
# Read your shapefile:
# Load the shapefile for ITL2 regions (replace 'path_to_shapefile' with the actual file path)
# https://ckan.publishing.service.gov.uk/dataset/international-territorial-level-2-january-2025-boundaries-uk-bgc/resource/78abd24e-5f15-4e04-a118-ddc649e9c769
itl2_shapefile <- st_read("~/mSCAPE/1_epidemic_modelling/UK_region_lookups/ITL2 shape files/ITL2_JAN_2025_UK_BFE.shp")

# Convert migration matrix to a long format data frame:
#migration_df <- as.data.frame( commuting_ITL2_inprob_list[[1]] )
commuting_ITL2_inprob_df$origin <- rownames( commuting_ITL2_inprob_df ) #** Check which way around for origin and destination **
commute_long <- pivot_longer( commuting_ITL2_inprob_df,
                               cols = -origin,
                               names_to = "destination",
                               values_to = "flow")
# Join the migration data with your shapefile:
uk_map_with_flows <- left_join( itl2_shapefile, commute_long,
                                  by = c("ITL225CD" = "origin"))
# Create the static map:
commute_flow_plot <- ggplot(data = uk_map_with_flows) +
                          geom_sf(aes(fill = flow)) +
                          scale_fill_viridis_c(option = "plasma", name = "Commuting Flow") +
                          theme_minimal() +
                          labs(title = "Commuting Flows Map",
                               subtitle = "Based on ONS UK 2021 census origin-destination data")
