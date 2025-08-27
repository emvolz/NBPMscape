#### Site selection using rational 'greedy' algorithm

# Author: Kieran Drake
# Date: 26 august 2025

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
region_TD <- read.csv("~/GitHub/NBPMscape/scripts/median_TD_by_region/sim_regionentry/60000_replicates/psampled_100/TD_by_region_sorted_simreps_60000_psampled_100.csv")

# Load correlations between region pair TD (time to detection) 
region_TD_corr <- read.csv("~/GitHub/NBPMscape/scripts/median_TD_by_region/sim_regionentry/60000_replicates/psampled_100/regional_TD_correlation_matrix_simreps_60000_psampled_100.csv")

corr_matrix = region_TD_corr[,c(2:40)]
colnames(corr_matrix) = region_TD_corr[,1]

# Define function for allocating sites to region while considering the block 
# they are in
site_select_by_block <- function(  region_td_data = region_TD
                                 , corr_mat = corr_matrix
                                 , n_blocks = 10
                                 , n_sites = 25
                                 )
{

  # Perform hierarchical clustering
  hc <- hclust( as.dist( 1 - abs(corr_mat) ), method = "average")
  # Cluster/block assignment for each variable
  block_membership <- as.data.frame( cutree( hc, k = n_blocks) )
  ITL2_code_name <- as.data.frame(rownames(block_membership))
  ITL2_code <- substr(ITL2_code_name[,1],1,4)
  block_membership <- cbind(ITL2_code,ITL2_code_name, block_membership)
  rownames(block_membership) <- NULL
  colnames(block_membership) <- c("ITL2_code","ITL2_name_code","block_membership")
  
  # Merge cluster/block assignment with time to detection data
  region_td_block <- merge(region_td_data[,c(1,2,6)], block_membership[,c(1,3)]
                          , by = "ITL2_code", all = TRUE)
  
  # Sort by time to detection
  region_td_block_sorted <- region_td_block[ order( region_td_block[,"treport_median_fit"], decreasing = FALSE ) , ]
  
  # Add column to store site allocations
  region_td_block_sorted[ , ncol( region_td_block_sorted ) + 1 ] <- rep(0,nrow(region_td_block_sorted))
  colnames(region_td_block_sorted)[ ncol( region_td_block_sorted ) ] <- "sites_allocated"
    
  ## Greedy algorithm
  
  # Initiate df to track the number of sites allocated per block
  block_regions_sites <- data.frame( block = seq(1,n_blocks)
                                    , regions = rep(0,n_blocks)
                                    , number_of_sites = rep(0,n_blocks))
  for(i in 1:n_blocks){
    block_regions_sites$regions[i] <- nrow(subset(region_td_block_sorted, block_membership == i))
  }
  
  # Set limit for the number of sites per region
  region_site_limit = ceiling( n_sites / nrow( region_td_block_sorted ) )
  
  # Initialise block allocation limit - starting at 1
  block_allocation_limit = 1
  
  # Repeat until all sites are allocated
  n_sites_unallocated <- n_sites
  
  # Define not in
  `%ni%` <- Negate(`%in%`)
  
  while (n_sites_unallocated > 0) {
    
    # Loop through regions from shortest TD to longest
    for (i in 1:nrow( region_td_block_sorted )){
      
      if( n_sites_unallocated == 0 ) {break}
      
      ############################
      # TESTING
      #i<-39
      #region_td_block_sorted[ i, ]
      ############################
      
      ### Apply rules for site allocation
      
      # Obtain list of blocks in which there are regions without a site allocated
      blocks_w_empty_regions <- subset( block_regions_sites, regions > number_of_sites )
      
      # 1
      # Move to next row (region) if ...
      # all regions in this particular region's block already has a site allocated
      if(
          region_td_block_sorted$block_membership[ i ] %ni% blocks_w_empty_regions$block 
      ){ 
        #print( "TRUE" )
        next 
        } #else { print("FALSE")}
      
      # 2
      # Move to next row if ... 
      # block already has more sites allocated than other blocks
      block_regions_sites_temp_subset <- block_regions_sites[ region_td_block_sorted$block_membership[ i ], ]
      if( 
        block_regions_sites_temp_subset$number_of_sites >= block_allocation_limit 
      ){ 
        #print( "TRUE" )
        next 
      } #else { print("FALSE")}
      
      # 3
      # Also need to exclude blocks where the number of sites allocated is 
      # already equal to the number of regions
      # SAME AS FIRST
      if( 
        block_regions_sites_temp_subset$number_of_sites == block_regions_sites_temp_subset$regions
      ){ #print( "TRUE" )
        next 
      } #else { print("FALSE")}
      
      # 4
      # Move to next row if ... 
      # number of sites allocated to region is already equal to the limit for an individual region
      if( 
        region_td_block_sorted$sites_allocated[i] >= region_site_limit
        ){ #print( "TRUE" )
        next 
      } #else { print("FALSE")}
      
      # 5
      # Move to next row if ... 
      # some region(s) in the block have fewer sites allocated to them
      if( 
        sum( subset( region_td_block_sorted
                       , block_membership == region_td_block_sorted$block_membership[ i ] )$sites_allocated < region_td_block_sorted$sites_allocated[ i ] ) >= 1 
      ){ #print( "TRUE" )
        next 
      } #else { print("FALSE")}
      
      # 6
      # Move to next row if ... 
      # other region(s) (regardless of block) have fewer number of sites allocated to them
      if( 
        sum( region_td_block_sorted$sites_allocated[i] > region_td_block_sorted$sites_allocated ) >= 1 
      ){ #print( "TRUE" )
        next 
      } #else { print("FALSE")}
      
      ## Add to region and block site allocation numbers
          
      # Select region for site and add to region site total
      region_td_block_sorted$sites_allocated[i] <- region_td_block_sorted$sites_allocated[i] + 1
      
      # Also add to total number of sites per block
      block_regions_sites[region_td_block_sorted$block_membership[i],3] <- block_regions_sites[region_td_block_sorted$block_membership[i],3] + 1 
      
      # Once site is allocated, reduce the number remaining to be allocated
      n_sites_unallocated <- n_sites_unallocated - 1
      
      ## Adjust allocation limits   
      
      # The number of sites to allocate to a single block is limited until all 
      # blocks reach that limit and then the limit is increased
      # Only increase the block site allocation limit if...
      
      # If the number of regions in a block is the same as the number of sites in
      # a block then increase the block allocation limit
      #if( sum( block_regions_sites$regions == block_regions_sites$number_of_sites ) == nrow(block_regions_sites) ){
      #  block_allocation_limit <- block_allocation_limit + 1
      #}
      
      # Exclude blocks that have a filled the allocation for their regions
      block_site_limits <- block_regions_sites$regions * region_site_limit
      block_regions_sites_ex_full_alloc <- subset( block_regions_sites, number_of_sites != block_site_limits )
      
      # If all blocks have been allocated the same number of sites, then increase
      # the block allocation limit
      if( ( length( unique( block_regions_sites_ex_full_alloc$number_of_sites ) ) == 1 ) & # ...the number of sites allocated to each block is the same, and...
            unique( block_regions_sites_ex_full_alloc$number_of_sites )[1] > 0 ){ # ...the number of sites in each block is greater than 0
          
        block_allocation_limit <- block_allocation_limit + 1     # Increase block allocation limit
        
        print( paste0( "block_allocation_limit = ", block_allocation_limit ) )
        
        break # Start from beginning of for loop once block_allocation_limit increased - otherwise might allocation to regions with longer time to detection before shorter TD
        
      } else{
            
        block_allocation_limit <- block_allocation_limit # No change to block allocation limit
        print( paste0( "block_allocation_limit = ", block_allocation_limit ) )      
      
      } 

    }
    
  }
return( list( region_td_block_sorted, block_regions_sites ) )  
}

# Run greedy algorithm with 10 blocks
site_select_output_10 <- site_select_by_block(  region_td_data = region_TD, corr_mat = corr_matrix
                                           , n_blocks = 10, n_sites = 25 )
# Run greedy algorithm with 20 blocks
site_select_output_20 <- site_select_by_block(  region_td_data = region_TD, corr_mat = corr_matrix
                                                , n_blocks = 20, n_sites = 25 )
# Run greedy algorithm with 30 blocks
site_select_output_30 <- site_select_by_block(  region_td_data = region_TD, corr_mat = corr_matrix
                                                , n_blocks = 30, n_sites = 25 )

## Checks
View( site_select_output_10[[1]] )
View( site_select_output_10[[2]] )

View( site_select_output_20[[1]] )
View( site_select_output_20[[2]] )

View( site_select_output_30[[1]] )
View( site_select_output_30[[2]] )

sum( site_select_output_10[[1]]$sites_allocated ) == n_sites
sum( site_select_output_20[[1]]$sites_allocated ) == n_sites
sum( site_select_output_30[[1]]$sites_allocated ) == n_sites