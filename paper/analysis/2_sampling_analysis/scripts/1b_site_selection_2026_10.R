#### Site selection using rational 'greedy' algorithm

# Author: Kieran Drake
# Date: 26 August 2025. Updated March 2026 and October 2026

## Idea is to:
# allocate sites to the regions with the shortest time to detection 
# AND
# to distribute evenly amongst blocks (which are formed of regions that are 
# similar, i.e. determined using correlation of times to detection series from 
# multiple simulation replicates)
# ALTHOUGH (in this case)
# blocks with more regions will be allocated more sites
# Implementation:
# 1) sort regions by metric (e.g. median time to detection)
# 2) allocate sites to the region with the best metric (top of the list, shortest
#    time to detection),
#    BUT check rules first
# Rules to follow:
# 3) allocate equal number of sites to each block before allocating more sites 
#    to an individual block (e.g. all blocks must have 1 site before any block 
#    can have 2 sites)
# 4) allocate equal number of sites to each region before allocating more sites
#    to an individual region (e.g. all regions must have 1 site before any region 
#    can have 2 sites)
# 5) some blocks have more regions than others and so this impacts the calculation
#    for rules 3 & 4 
#    (i.e. if block 7 has 9 regions and 2 sites allocated already and block 3 has
#    only 1 region with 1 site allocated then we would not allocate more sites to
#    block 3 before allocating a site to each of the 9 regions in block 7. So block 7
#    can have 9 sites allocated, possibly 10, before we add a 2nd site to block 3's 
#    only region)

# Load data containing region and values for parameter to be used to determine site selection
region_TD <- read.csv("~/GitHub/NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_td_regional/TD_by_region_sorted_simreps_1000_psampled_100.csv")

# Load correlations between region pair TD (time to detection) 
region_TD_corr <- read.csv("~/GitHub/NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_td_regional/regional_TD_correlation_matrix_simreps_1000_psampled_100.csv")

corr_matrix = region_TD_corr[,c(2:40)]
#colnames(corr_matrix) = region_TD_corr[,1]

#** Remove ITL2 regions in Wales as these cannot currently be sampled**
wales_itl2 = c("TLL3", "TLL4", "TLL5")
# Find cols and rows corresponding to ITL2 regions in Wales
wales_itl2_indices_cm = c()
wales_itl2_indices_td = c()
for (i in 1:length(wales_itl2)){ #i=1
  wales_itl2_i = wales_itl2[i]
  wales_itl2_i_index_cm <- which( stringr::str_detect(string = names(corr_matrix)
                                                  , pattern = wales_itl2_i))
  wales_itl2_i_index_td <- which( stringr::str_detect(string = region_TD$ITL2_code
                                                      , pattern = wales_itl2_i))
  wales_itl2_indices_cm = append(wales_itl2_indices_cm, wales_itl2_i_index_cm)
  wales_itl2_indices_td = append(wales_itl2_indices_td, wales_itl2_i_index_td)
}
# Remove from correlation matrix
corr_matrix <- corr_matrix[-wales_itl2_indices_cm,-wales_itl2_indices_cm]
# Remove from TD df
region_TD <- region_TD[-wales_itl2_indices_td,]

# Define function for allocating sites to region while considering the block 
# they are in
#** Version 2 - updated to allow number of sites to be greater than the number of regions **
site_select_by_block_v2 <- function(  region_td_data = region_TD
                                      , corr_mat = corr_matrix
                                      , n_blocks = 10 #5
                                      , n_sites = 37 #40
)
{
  
  # Perform hierarchical clustering
  hc <- hclust( as.dist( 1 - abs(corr_mat) ), method = "average")
  # Cluster/block assignment for each variable
  block_membership <- as.data.frame( cutree( hc, k = n_blocks) )
  #ITL2_code_name <- as.data.frame(rownames(block_membership))
  ITL2_code_name <- as.data.frame(colnames(corr_mat))
  ITL2_code <- substr(ITL2_code_name[,1],1,4)
  block_membership <- cbind(ITL2_code,ITL2_code_name, block_membership)
  rownames(block_membership) <- NULL
  colnames(block_membership) <- c("ITL2_code","ITL2_name_code","block_membership")
  
  # Merge cluster/block assignment with time to detection data
  region_td_block <- merge(region_td_data[,c(1,2,6)], block_membership[,c(1,3)]
                           , by = "ITL2_code", all = TRUE)
  
  # Sort by time to detection
  region_td_block_sorted <- region_td_block[ order( region_td_block[,"treport_median_fit"], decreasing = FALSE ) , ]
  rownames(region_td_block_sorted) <- NULL
  
  # Add column to store site allocations
  region_td_block_sorted[ , ncol( region_td_block_sorted ) + 1 ] <- rep(0,nrow(region_td_block_sorted))
  colnames(region_td_block_sorted)[ ncol( region_td_block_sorted ) ] <- "sites_allocated_to_region"
  
  ## Greedy algorithm
  
  # Initiate df to track the number of sites allocated per block
  block_regions_sites <- data.frame( block = seq(1,n_blocks)
                                     , regions = rep(0,n_blocks)
                                     , number_of_sites = rep(0,n_blocks))
  for(i in 1:n_blocks){
    block_regions_sites$regions[i] <- nrow(subset(region_td_block_sorted, block_membership == i))
  }
  
  # Set limit for the number of sites per region
  region_site_limit_max = ceiling( n_sites / nrow( region_td_block_sorted ) )
  region_site_limit_running = 1 # Start at 1 and increase by 1 when all regions have the same number of sites, with the maximum being 'region_site_limit_max'
  
  # Initialise block allocation limit - starting at 1
  block_allocation_limits = sort(unique(block_regions_sites$regions))
  block_allocation_limit_running = block_allocation_limits[1]
  
  # Repeat until all sites are allocated
  n_sites_unallocated <- n_sites
  
  # Set the total number of sites allocatd to the block to 0 to begin with
  region_td_block_sorted$sites_allocated_to_block = 0
  
  # Add columns to df to track site allocations and limits
  region_td_block_sorted$current_region_site_limit <- rep(1,nrow(region_td_block_sorted))
  
  for (i in 1:nrow(region_td_block_sorted)){ #i=1
    region_td_block_sorted$current_block_site_limit[i] <- nrow(subset(region_td_block_sorted, block_membership==block_membership[i]))
  }

  block_regions_sites$current_block_sites_limit = block_regions_sites$regions * region_site_limit_running
  
    # Define not in
  `%ni%` <- Negate(`%in%`)
  
  while (n_sites_unallocated > 0) {
    
    for (i in 1:nrow(region_td_block_sorted)){ # i=1
      
      ### Assess region i for site allocation
      # Check criteria for site allocation
      # If any criteria are not met then move to next region (row) for assessment of site allocation
      
      ## Define values to evaluate criteria for region i and the block it is in
      sites_allocated_to_region_i <- region_td_block_sorted$sites_allocated_to_region[i]      
      current_region_site_limit_i <- region_td_block_sorted$current_region_site_limit[i]
      
      total_sites_allocated_to_region_i_block <- subset(block_regions_sites
                                                  , block == region_td_block_sorted$block_membership[i] )$number_of_sites
      current_sites_limit_region_i_block <- subset(block_regions_sites
                                                  , block == region_td_block_sorted$block_membership[i] )$current_block_sites_limit
      
      sites_allocated_across_region_i_block_df <- subset( region_td_block_sorted
                                                    , block_membership==block_membership[i])$sites_allocated_to_region
      region_i_current_block_site_limit <- region_td_block_sorted$current_block_site_limit[i]
      
      # Blocks with capacity for more sites to be allocated within current block site limits      
      blocks_w_capacity <- subset( block_regions_sites
                                   ,(block_regions_sites$current_block_sites_limit - block_regions_sites$number_of_sites) > 0)
      
      
      if( sites_allocated_to_region_i >= current_region_site_limit_i ){ #region_site_limit_running ){next}
        # Already at current region site limit
        next
      } else if( total_sites_allocated_to_region_i_block >= current_sites_limit_region_i_block ){
        #sum( subset( region_td_block_sorted, block_membership==block_membership[i])$sites_allocated_to_region ) >= block_allocation_limit_running ){
        # Already at site allocation limit for block
        next
      } else if( sites_allocated_to_region_i > min( sites_allocated_across_region_i_block_df ) ){
        # Would not take region more than 1 site above other regions in same block
        next
      } else if(
          # Would not take region i's block total sites to more than 1 above blocks 
          # that still have capacity to add sites in this round (i.e. ensuring that 
          # each region, regardless of block, has the same number of sites before 
          # additional sites are allocated to a block)
          # 
          #(region_td_block_sorted$current_block_site_limit[i] - region_td_block_sorted$sites_allocated_to_block[i]) >
          #          min(blocks_w_capacity$current_block_sites_limit - blocks_w_capacity$number_of_sites)
          #(region_i_current_block_site_limit - total_sites_allocated_to_region_i_block) >
          #    min(blocks_w_capacity$current_block_sites_limit - blocks_w_capacity$number_of_sites)
          
          # If the difference between the number of sites allocated to any block and 
          # the number allocated to the block containing region i plus 1 is greater 
          # than 1, then move to next region in list sorted by shortest TD
          # Number of sites allocatd to region i's block + 1 > 
          sum( ((total_sites_allocated_to_region_i_block+1) - blocks_w_capacity$number_of_sites ) > 1) > 0 ){
        next
        }
      
      ### All criteria met for site allocation, so add a site
      region_td_block_sorted$sites_allocated_to_region[i] <- region_td_block_sorted$sites_allocated_to_region[i] + 1
      print(region_td_block_sorted[i,])
      
      # Sites allocated to block
      regions_in_block <- subset(region_td_block_sorted, block_membership == region_td_block_sorted$block_membership[i] )
      sites_in_block <- sum(regions_in_block$sites_allocated_to_region)
      for ( j in 1:nrow(region_td_block_sorted)){
        if(region_td_block_sorted$block_membership[j] == region_td_block_sorted$block_membership[i]){
          region_td_block_sorted$sites_allocated_to_block[j] = sites_in_block  
        }  
      }
      #block_regions_sites$number_of_sites = c(1,0,1,1,0)
      for (k in 1:nrow(block_regions_sites)){ # k=3
        block_n = block_regions_sites$block[k]
        block_regions_sites$number_of_sites[k] = sum( subset( region_td_block_sorted
                                                           , block_membership == block_n)$sites_allocated_to_region )
      }
      
      ### Re-assess limits on block, region and site allocation limits
      
      # Reduce number of sites left to allocate
      n_sites_unallocated <- n_sites_unallocated - 1
      # If all sites have been allocated then exit loop
      if( n_sites_unallocated == 0 ){ break }
      
      # Increase region site limit if all regions at the current limit of sites allocated
      if( nrow( subset( region_td_block_sorted, sites_allocated_to_region == region_site_limit_running )) == nrow(region_td_block_sorted) ){
        region_site_limit_running <- min( region_site_limit_running + 1, region_site_limit_max) # But do not increase above maximum limit
      
        # Update limits in main df
        # Region site limits
        region_td_block_sorted$current_region_site_limit = region_site_limit_running

        # Block site limits
        for (i in 1:nrow(region_td_block_sorted)){ #i=1
          region_td_block_sorted$current_block_site_limit[i] <- nrow(subset(region_td_block_sorted, block_membership==block_membership[i])) * region_site_limit_running
        }
        block_regions_sites$current_block_sites_limit = block_regions_sites$regions * region_site_limit_running
      
        # Blocks with capacity for more sites to be allocated within current block site limits      
        blocks_w_capacity <- subset( block_regions_sites
                                     ,(block_regions_sites$current_block_sites_limit - block_regions_sites$number_of_sites) > 0)
        
        #... and restart at first row (region with shortest TD)
        break
      }
      
      # Blocks with capacity for more sites to be allocated within current block site limits      
      blocks_w_capacity <- subset( block_regions_sites
                                   ,(block_regions_sites$current_block_sites_limit - block_regions_sites$number_of_sites) > 0)
      
      # Restart for loop every time allocate a site
      break
      
      #if(nrow(blocks_w_capacity) == 1)
        
      # If all blocks have 
      
      # if all blocks with 1 region have 1 site then increase block_allocation_limit to 2
      # then, if all blocks with 2 regions have 2 sites (i.e. one in each region) then increase block_allocation_limit to next level (in 5block model this would be 31 and
      # would then allow allocation of a site to each region in 31-region block)
      #block_regions_sites
      #if( subset(region_td_block_sorted$))
        
      
      
    }
    print(n_sites_unallocated)
    
  }
  return( list( region_td_block_sorted, block_regions_sites ) )  
}

#### Run greedy algorithm to allocate sites using different numbers of blocks
#### and different numbers of sites
#### (this section was updated using Claude Sonnet 5 via Imperial College dAIsy on 16 Oct 2026)

# ---- Define combinations of blocks and sites ----
param_grid <- expand.grid(
  n_blocks = c(5, 10, 20, 30),
  n_sites  = c(8, 10, 20, 26, 29, 30, 40, 50, 60),
  KEEP.OUT.ATTRS   = FALSE,
  stringsAsFactors = FALSE
)

# ---- Run 'site_select_by_block_v2' function for all block/sites combinations (once per row of param_grid) ----
site_select_outputs <- Map(
  function(n_blocks, n_sites) {
    site_select_by_block_v2(
      region_td_data = region_TD,
      corr_mat       = corr_matrix,
      n_blocks       = n_blocks,
      n_sites        = n_sites
    )
  },
  param_grid$n_blocks,
  param_grid$n_sites
)

# ---- Name results with n block and m site labels ----
names(site_select_outputs) <- paste0(
  "site_select_output_", param_grid$n_blocks, "blocks_",
  param_grid$n_sites, "sites"
)

# ---- Access results if required, e.g.: ----
# site_select_outputs$site_select_output_5blocks_8sites_v2
# site_select_outputs$site_select_output_10blocks_10sites_v2
# site_select_outputs$site_select_output_20blocks_20sites_v2

# ---- (Optional) unpack into individual objects in your environment ----
# list2env(site_select_outputs, envir = .GlobalEnv)

# ---- Create an output directory to store the CSV files ----
output_dir <- file.path(normalizePath("~", mustWork = FALSE), "GitHub/NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/site_select_outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Add sites_allocated column (copied from sites_allocated_to_region)
# ---- Check total allocated sites matches n_sites and print summary at end
# ---- Write each output's first element to a .csv file named after the output

# Create object to hold test results for each combination of blocks and sites
test_results <- vector("list", length(site_select_outputs))

# Loop through all outputs (combinations of blocks and sites)
for (i in seq_along(site_select_outputs)) {#i=1
  output_name <- names(site_select_outputs)[i]
  n_sites_i   <- param_grid$n_sites[i]
  
  # Add 'sites_allocated' column
  site_select_outputs[[output_name]][[1]]$sites_allocated <-
    site_select_outputs[[output_name]][[1]]$sites_allocated_to_region
  
  # ---- Check: total sites allocated to region should equal n_sites ----
  actual_sum   <- sum(site_select_outputs[[output_name]][[1]]$sites_allocated_to_region)
  check_passed <- actual_sum == n_sites_i
  
  test_results[[i]] <- data.frame(
    output_name = output_name,
    n_blocks    = param_grid$n_blocks[i],
    n_sites     = n_sites_i,
    actual_sum  = actual_sum,
    passed      = check_passed,
    stringsAsFactors = FALSE
  )
  
  if (check_passed) {
    write.csv(
      site_select_outputs[[output_name]][[1]],
      file      = file.path(output_dir, paste0(output_name, ".csv")),
      row.names = FALSE
    )
  }
}

# ---- Combine into a single summary data frame and print once at the end ----
test_summary <- do.call(rbind, test_results)
print(test_summary)

# ---- Flag loudly if any check failed (number of sites allocated does not match intended number of site),
#      so it can't be missed ----
if (!all(test_summary$passed)) {
  warning(
    "One or more checks failed — see test_summary for details:\n",
    paste(
      sprintf(
        "  %s (n_sites expected = %d, actual_sum = %d)",
        test_summary$output_name[!test_summary$passed],
        test_summary$n_sites[!test_summary$passed],
        test_summary$actual_sum[!test_summary$passed]
      ),
      collapse = "\n"
    )
  )
}

# ---- If want to view any of the site allocation results
#View(site_select_outputs$site_select_output_5blocks_8sites[[1]])
#View(site_select_outputs$site_select_output_5blocks_8sites[[2]])


### Create df containing ITL2 code, ITL2 name, treport_median_fit, sites/blocks as
# columns and populated with the number of sites allocated

# Create base
site_select_df <- site_select_output_5blocks_8sites_v2[[1]]
# Trim to desired columns
site_select_df <- site_select_df[,c("ITL2_code","ITL2_name","treport_median_fit")]

# Loop through individual dataframes and copy column containing number of sites 
# allocated to the compilation df
for (i in seq_along(site_select_outputs)) {#i=2
  output_name <- names(site_select_outputs)[i]
  
  # Add column from df for specific block/site combination to compilation df
  site_select_df$new_col <- site_select_outputs[[output_name]][[1]]$sites_allocated_to_region #temp_df$sites_allocated_to_region
  
  # Rename new column added
  # Extract number of blocks and sites from the output_name
  n_blocks <- regmatches( output_name, regexpr("(?<=_)(\\d+)(?=blocks)", output_name, perl = TRUE) )
  n_sites  <- regmatches( output_name, regexpr("(?<=_)(\\d+)(?=sites)",  output_name, perl = TRUE) )
  
  new_name <- paste0( n_blocks,"b_",n_sites,"s")
  names(site_select_df)[names(site_select_df) == "new_col"] <- new_name
  
}

# Write to file
write.csv( site_select_df, file = file.path( output_dir, paste0("site_select_output_compilation.csv") ), row.names = FALSE)
