#=
Function for selecting sampling sites

=#


"""
Function    select_nhs_trust_sampling_sites

Description     Function to select sites (NHS Trusts) that will be used for sampling.
                See 'Examples' below for format of input file and output file.

Arguments   prioritise_adult_sampling::Boolean      Function selects sites based on block allocation, and probability of
                                                    attending a particular NHS Trust given the ITL2 region. These probabilities
                                                    are disaggregated for adults and children and so the function alternates
                                                    between adult or child probabilities to rank NHS Trusts for selection.
                                                    If 'true' then the first site will be based on adult ranking, if 'false'
                                                    then it will be based on child ranking.
            n_sites_max::Int                        Maximum number of sites (127 is the number of NHS Trusts in England relevant for ARI)
            site_selection_folder::String           Folder containing the number of sites allocated to each block
            site_selection_file::String             Name of file containing the number of sites allocated to each block
            output_folder::String                   Folder where output file should be saved
            output_file::String                     Output file name

Returns     Writes a dataframe to the file and folder specified.
            Output file will contain columns:
                    NHS_Trust_code, ICU_site_name, ICU_full_name,
                    ICU_beds_adult,	ICU_beds_paediatric,
                    sample_target_per_year, sample_target_per_week

Examples    
            # Input file of format below, as output from 'site_selection_2026_03.R' 
            # (first three rows from site selection through hierarchical clustering and geographic stratification (HCGS)
            # using 20 blocks and 26 sites)
            ITL2_code	ITL2_name	                        treport_median_fit	block_membership	sites_allocated_to_region	sites_allocated_to_block	current_region_site_limit	current_block_site_limit	sites_allocated
            TLI7	    Outer London - West and North West	45.79655696	        18	                1	                        1	                        1	                        1	                        1
            TLJ2	    Surrey, East and West Sussex	    46.58666992	        4	                1	                        2	                        1	                        4	                        1
            TLD3	    Greater Manchester	                47.15547918	        4	                1	                        2	                        1	                        4	                        1
            ...         ...                                 ...                 ...                 ...                         ...                         ...                         ...                         ...
            
            # Run function
            select_nhs_trust_sampling_sites(; prioritise_adult_sampling = true
                                                , n_sites_max = 127
                                                , site_selection_folder = joinpath( @__DIR__, "..", "examples/hcgs_site_selection")
                                                , site_selection_file = "example_site_select_output_20blocks_26sites.csv"
                                                , output_folder = joinpath( @__DIR__, "..", "examples/hcgs_site_selection")
                                                , output_file = "example_sampling_sites_20blocks_26sites.csv"
                                                )
            
            

"""
function select_nhs_trust_sampling_sites(; prioritise_adult_sampling = true
                                            , n_sites_max = 127
                                            , site_selection_folder
                                            , site_selection_file
                                            , output_folder
                                            , output_file
                                            )
    # Return list of NHS Trusts using the number of sites allocated to ITL2 regions
    n_sites_by_region = CSV.read( joinpath( site_selection_folder, site_selection_file), DataFrame ) #site_select_output_10blocks_10sites

    #prioritise_adult_sampling = true
    #n_sites_max = 127 # The total number of NHS Trusts

    nhs_trust_site_list = []#Vector{String}
    for i in 1:nrow(n_sites_by_region) #i=34
        itl2_code = n_sites_by_region[i,:ITL2_code]
        n_sites = n_sites_by_region[i,:sites_allocated]
        n_sites_allocated = 0
        allocation_attempts = 0
        if n_sites > 0
            nhs_trust_probs_adult = sort(ITL2_TO_NHS_TRUST_PROB_ADULT[:,[:NHS_Trust_code,Symbol(itl2_code)]], Symbol(itl2_code), rev=true)
            nhs_trust_probs_child = sort(ITL2_TO_NHS_TRUST_PROB_CHILD[:,[:NHS_Trust_code,Symbol(itl2_code)]], Symbol(itl2_code), rev=true)

            # Helper function to determine the order in which sites are allocated to NHS Trusts based on adult or paediatric catchments.
            # If adult NHS Trust allocation is prioritised then odd number allocation attempts will be adult NHS Trusts and paediatric will be even and vice versa
            function adult_child_site_selection_order(; allocation_attempts = allocation_attempts, prioritise_adult_sampling = prioritise_adult_sampling)
                if prioritise_adult_sampling == true
                    result = isodd( allocation_attempts )
                elseif prioritise_adult_sampling == false
                    result = iseven( allocation_attempts )
                end
                return( result ) 
            end
                        
            while n_sites_allocated != n_sites # 
                
                allocation_attempts = allocation_attempts +1

                if adult_child_site_selection_order(; allocation_attempts = allocation_attempts, prioritise_adult_sampling = prioritise_adult_sampling) 

                    # Check if top ranked NHS Trust by probability is already in the list
                    # and if it is then remove it from the list of possible NHS Trusts to add to list
                    while nhs_trust_probs_adult[1,:NHS_Trust_code] in nhs_trust_site_list
                        nhs_trust_probs_adult = nhs_trust_probs_adult[2:end,:]
                    end
                    # Check that NHS Trust has non-zero probability of attendance
                    # If it does then add it to the list...
                    if nhs_trust_probs_adult[1,2] > 0.0
                        push!( nhs_trust_site_list, nhs_trust_probs_adult[1,:NHS_Trust_code])
                        n_sites_allocated = n_sites_allocated + 1
                    end

                else

                    while nhs_trust_probs_child[1,:NHS_Trust_code] in nhs_trust_site_list#println(nhs_trust_site_list)
                        nhs_trust_probs_child = nhs_trust_probs_child[2:end,:]
                    end
                    
                    # Check that NHS Trust has non-zero probability of attendance
                    # If it does then add it to the list...
                    if nhs_trust_probs_child[1,2] > 0.0
                        push!( nhs_trust_site_list, nhs_trust_probs_child[1,:NHS_Trust_code])
                        n_sites_allocated = n_sites_allocated + 1
                    end
                end
            end
        end
    end


    # Generate file of selected sites (NHS Trusts) and their adult and paediatric ICU beds
    sampling_sites_list = filter( row -> row.NHS_Trust_code in nhs_trust_site_list, ARI_CC_BED_SITREP[:, [:NHS_Trust_code
                                                                                                        ,:NHS_Trust_name
                                                                                                        ,:Adult_critical_care_beds_available
                                                                                                        ,:Paediatric_intensive_care_beds_available] ] )
    # Adjust so matches the format required for icu_td function
    # Target column names = NHS_Trust_code, ICU_site_name, ICU_full_name, ICU_beds_adult, ICU_beds_paediatric, sample_target_per_year,sample_target_per_week
    rename!(sampling_sites_list, :NHS_Trust_name => :ICU_site_name
                                , :Adult_critical_care_beds_available => :ICU_beds_adult
                                , :Paediatric_intensive_care_beds_available => :ICU_beds_paediatric
    )
    insertcols!(sampling_sites_list, 3, :ICU_full_name => Vector{Union{Missing, String}}(missing, nrow(sampling_sites_list)))
    insertcols!(sampling_sites_list, 6
                ,:sample_target_per_year => Vector{Union{Missing, Int}}(missing, nrow(sampling_sites_list))
                ,:sample_target_per_week => Vector{Union{Missing, String}}(missing, nrow(sampling_sites_list))
                )

    # Check that sampling sites selected total the number required
    if nrow(sampling_sites_list) == sum(n_sites_by_region.sites_allocated)
        println( "Sites required: $(sum(n_sites_by_region.sites_allocated)), sites selected: $(nrow(sampling_sites_list))")
    else
        @warn "The number of sites required is different from the number selected! Sites required: $(sum(n_sites_by_region.sites_allocated)), sites selected: $(nrow(sampling_sites_list))"
    end

    CSV.write( joinpath( output_folder, output_file), sampling_sites_list )
end