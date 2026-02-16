File: script.R
Description:	Inspecting the origin-destination data from the UK 2021 census,
		published by the Office for National Statistics at https://www.nomisweb.co.uk/sources/census_2021_od
		Script uses output from region_mapping.R and outputs commuting_ITL2_matrix.rds.

File: region_mapping.R
Description: Mapping from UTLA to ITL2

File: prepmig.R
Description: Formats migration data for inclusion in Julia code. Converts commuting_ITL2_matrix.rds to:
		- commuting_ITL2_list.rds
		- commuting_ITL2_inprob_list.rds
		- commuting_ITL2_prob_list.rds

File: ITL2_key2.rds
Decsription: ITL2 region codes and names in England and Wales


Data files used in 'region_mapping.R':

###########################################################
# ITL2
ITL2_JAN_2025_UK_BFC_-5615416600314916039.csv
# Downloaded from ONS on 17 March 2025
https://geoportal.statistics.gov.uk/datasets/f2f482dc109e454599829bc5ab2c4418_0/explore
https://geoportal.statistics.gov.uk/datasets/ons::international-territorial-level-2-january-2025-boundaries-uk-bfc/about

###########################################################
Data set for mapping from LTLA to ITL2
LAD_(December_2024)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2025)_Lookup_in_the_UK.csv
Downloaded from ONS on 17 March 2025
https://geoportal.statistics.gov.uk/datasets/15bdb6b0ff1c4a64b34a64e2e39f8caf_0/explore

############################################################
Mapping changes for ITL2 between 2021 and 2025
ITL_Level_2_(2021)_to_ITL_Level_2_(2025)_Lookup_in_the_UK.csv
Downloaded from ONS on 17 March 2025
https://geoportal.statistics.gov.uk/datasets/ons::itl-level-2-2021-to-itl-level-2-2025-lookup-in-the-uk/about

#############################################################
Areas as at 2020 and 2021
LAD20_LAU121_ITL321_ITL221_ITL121_UK_LU_v2.xlsx
Downloaded from ONS on 17 March 2025
https://geoportal.statistics.gov.uk/documents/ons::lad-2020-to-lau1-to-itl3-to-itl2-to-itl1-january-2021-lookup-in-the-uk-v2/about

#############################################################
Lookup for MSOA changes between 2011 and 2021 (the times of the ONS censuses)
MSOA_(2011)_to_MSOA_(2021)_to_Local_Authority_District_(2022)_Best_Fit_Lookup_for_EW_(V2).csv
Downloaded from ONS on 8 October 2025
https://www.data.gov.uk/dataset/9aa055b1-9b55-4d40-9be0-70bab8d55889/msoa-2011-to-msoa-2021-to-local-authority-district-2022-best-fit-lookup-for-ew-v21

###########################################################
Mapping for various UK geo areas
pcd_oa_lsoa_msoa_ltla_utla_rgn_ctry_ew_may_2021_lu.zip

Downloaded from ONS on 15 March 2025
https://open-geography-portalx-ons.hub.arcgis.com/datasets/ons::postcode-may-2021-to-oa-2021-to-lsoa-to-msoa-to-ltla-to-utla-to-rgn-to-ctry-may-2021-best-fit-lookup-in-ew/about

Postcode (May 2021) to OA (2021) to LSOA to MSOA to LTLA to UTLA to RGN to CTRY (May 2021) Best Fit Lookup in EW

Summary
Postcode Lookup

A best-fit lookup between postcodes, 2021 Census Output Areas (OA), LSOAs, MSOAs, LTLAs, UTLAs, regions (RGN), countries (CTRY) and nationalities (NAT) as at May 2021 in England and Wales. Postcodes are best-fitted by plotting the location of the postcode's mean address into the areas of the output geographies. (File size 23 MB).
Field Names - PCD, GRIDGB1E, GRIDGB1N, OA21CD, LSOA21CD, LSOA21NM, MSOA21CD, MSOA21NM, LTLA22CD, LTLA22NM, LTLA22NMW, UTLA22CD, UTLA22NM, UTLA22NMW, RGN22CD, RGN22NM, RGN22NMW, CTRY22CD, CTRY22NM, CTRY22NMW, NAT22CD, NAT22NM, NAT22NMW
Field Types - All Text
Field Lengths - 7, 6, 7, 9, 9, 60, 9, 60, 9, 60, 60, 9, 60, 60, 9, 30, 30, 9, 10, 10, 9, 20, 20

Details
Data File
CSV Collection
October 9, 2017
Info Updated
December 18, 2024
Data Updated
August 11, 2022 at 12:00 PM
Published Date
Public
Anyone can see this content
Custom License
View license details
Relevant Area

################################################

