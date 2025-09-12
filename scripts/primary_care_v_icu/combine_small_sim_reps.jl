#= This file contains three functions:
1   {combine_sim_reps}
    This function can be used to read in multiple simulation files and combine them for analysis of a larger group of
    simulation replicates

2(a) {n_sims_w_results}
     Takes one or more simulation files as input and returns the percentage of simulation replicates in those files that have
     at least 1 non-infinite value for the range of variables selected
2(b) {n_sims_w_results_icu_gp}
     Takes only two files as input - one pre-filtered for GP cases and one pre-filtered for ICU cases

and also
3   Code to try if having issues with loading .jld2 files, try looking at different package versions
=#


using JLD2
using DataFrames

#=1   {combine_sim_reps}
    This function can be used to read in multiple simulation files and combine them for analysis of a larger group of
    simulation replicates =#
function combine_sim_reps(;  sim_input_folders # in vector format where multiple paths
                            , sim_object_name = "sims" # "sims_G_icu_filter" # "sims_G_gp_filter"
                            , nrep = 1000
                            #, combined_sim_output_folder
                            #, combined_sim_output_filename
                          )
    sim_files = []
    for i in 1:length(sim_input_folders)
        sim_files_temp = readdir(sim_input_folders[i]; join=true)
        sim_files = vcat(sim_files, sim_files_temp)
    end

    # Create vector to store filtered dataframes
    sims = [DataFrame() for _ in 1:length(nrep)]

    # Add all simulation replicates from multiple files to a single object 
    for j in 1:length(sim_files)
        println("Processing ", sim_files[j])
        # Load simulation replicates from file
        sims_temp = load(sim_files[j], sim_object_name)

        # Number of replicates 
        nrep = length(sims_temp)
    
        # Append sim replicates from each file to a single vector with a df for each sim replicate
        if j == 1
            sims = sims_temp
        else
            sims = append!(sims, sims_temp)
        end

    end

    # Check length of file in terms of sim reps now
    if  nrep * length(sim_files) == length(sims)
        println("All simulation replicates added")
    end
    
    return(sims)

end

## Combine and save multiple files each containing 1000 simulation replicates that have been filtered to retain ICU cases only (G dataframe)
sims_G_icu_filter = combine_sim_reps( sim_input_folders = [ "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_icu"]
                                       #                     ,"C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_gp"]
                                    , sim_object_name = "sims_G_icu_filter" #"sims_G_gp_filter"# "sims" #  
                                    , nrep = 1000
                                    )
@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2" sims_G_icu_filter #@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter)_955898.jld2" sims_G_icu_filter

## Combine and save multiple files each containing 1000 simulation replicates that have been filtered to retain GP cases only (G dataframe)
sims_G_gp_filter = combine_sim_reps( sim_input_folders = [ "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_gp"]
                                    , sim_object_name = "sims_G_gp_filter" #"sims_G_icu_filter" ## "sims" #  
                                    , nrep = 1000
                                    )
@save "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2" sims_G_gp_filter #@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_gp_filter)_955898.jld2" sims_G_gp_filter

### Check whether can reload files - large files can be difficult to reload
sims_G_icu_filter = load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2", "sims_G_icu_filter") #load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2", "sims_G_icu_filter")
sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2", "sims_G_gp_filter") #sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2", "sims_G_gp_filter")

############################################
#### Count how many simulations contain results

# 2(a) {n_sims_w_results}
#      Takes one or more simulation files as input and returns the percentage of simulation replicates in those files that have
#      at least 1 non-infinite value for the range of variables selected
function n_sims_w_results(; sims_files = [], sims_object_names = [], count_types = [] ) # Input vector of files
    
    # Check whether the lengths match
    if length(sims_files) != length(sims_object_names)
        "The number of files and object names must be the same"
        return
    end

    # Create df to hold results
    results_simrep_perc_df = DataFrame()
    column_names = vcat("sims_files",count_types)
    num_rows = length(sims_files)
    columns = Dict(
        column_names[1] => Vector{String}(undef, num_rows),
        [column_names[i] => Vector{Float64}(fill(Inf, num_rows)) for i in 2:length(column_names)]...
    )
    results_simrep_perc_df = DataFrame(columns)

    # Loop through sim files
    for sf in 1:length(sims_files)
        
        # Load file
        sims_temp = load(sims_files[sf], sims_object_names[sf])
        nrep_sims = length(sims_temp)    
        
        # Initialise counters
        counter = zeros(Int, length(count_types))
        
        # Loop through simulation replicates in file 
        for sr in 1:length(sims_temp)
            
            
            # Check whether sim df is empty
            if size(sims_temp[sr],1) > 0
                # Loop through object types to count simulation replicates that have results for each count_type
                for c in 1:length(count_types)
                    # Filter sim df for finite values
                    sims_filtered =  sims_temp[sr][ isfinite.(sims_temp[sr][:,count_types[c]]), : ]
                    # Add to counter
                    if size( sims_filtered,1) > 0
                        counter[c] = counter[c] +1
                    end
                end
            end
        end
        
        #println("For sim file $(sims_files[sf]), the proportion of simulation replicates that have time values are:")
        
        #for ct in 1:length(count_types)
        #    println("   $(count_types[ct]): $( 100* counter[ct] / nrep_sims )%")
        #end
        results_simrep_perc_df[ sf , :] = vcat(sims_files[sf], counter / nrep_sims)        
    end

    # Print results
    println(results_simrep_perc_df)

end

n_sims_w_results(; sims_files = ["covidlike-1.3.1-sims-nrep1000_filtered_G_icu_955898.2.jld2"
                                ,"covidlike-1.3.1-sims-nrep1000_filtered_G_icu_955898.3.jld2"] #[]
                 , sims_object_names = [ "sims_G_icu_filter"
                                        ,"sims_G_icu_filter"] #[]
                 , count_types = [ "tgp", "ticu" ]
                 ) # Input vector of files

##Row │ sims_files                                                     tgp      ticu    
#     │ String                                                         Float64  Float64
#─────┼─────────────────────────────────────────────────────────────────────────────────
#   1 │ covidlike-1.3.1-sims-nrep1000_filtered_G_icu_955898.2.jld2    0.717    0.795
#   2 │ covidlike-1.3.1-sims-nrep1000_filtered_G_icu_955898.2.jld2    0.704    0.78


# 2(b) {n_sims_w_results_icu_gp}
#      Takes only two files as input - one pre-filtered for GP cases and one pre-filtered for ICU cases

function n_sims_w_results_icu_gp(; sims_file_icu, sims_file_gp)
    # Load files   
    sims_G_icu_filter = load(sims_file_icu, "sims_G_icu_filter")
    sims_G_gp_filter = load(sims_file_gp, "sims_G_gp_filter")
    
    # Count number of simulations
    nrep_sims_G_gp_filter = length(sims_G_gp_filter)
    nrep_sims_G_icu_filter = length(sims_G_icu_filter)

    # Initialise counters
    gp_results_counter = 0; icu_results_counter = 0

    # Count simulations with GP results
    for s in 1:nrep_sims_G_gp_filter
        if size(sims_G_gp_filter[s],1) >0
            gp_results_counter = gp_results_counter +1
        end
    end
    println("For GP cases, $(gp_results_counter) ($(round(100*gp_results_counter/nrep_sims_G_gp_filter))%) out of $(nrep_sims_G_gp_filter) simulations contain results")
    # Count simulations with ICU results
    for s in 1:nrep_sims_G_icu_filter
        if size(sims_G_icu_filter[s],1) >0
            icu_results_counter = icu_results_counter +1
        end
    end
    println("For ICU cases, $(icu_results_counter) ($(round(100*icu_results_counter/nrep_sims_G_icu_filter))%) out of $(nrep_sims_G_icu_filter) simulations contain results")
end

n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep53000.jld2" 
                , sims_file_gp = "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep53000.jld2") 
# ICU: 90% and GP: 90% 

n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep84000.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep84000.jld2") 
# Combination of HPC run 851401+854085 
# ICU: 35% and GP: 34% (because R was too low at ~1)

n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2") 
# ICU: 80% and GP: 80%


##########################################
# 3 Code to try if having issues with loading .jld2 files, try looking at different package versions

using Pkg
Pkg.status("JLD2")
println(Pkg.installed()["JLD2"])

using JLD2
Pkg.dependencies()[Base.UUID("e6f89c97-0c3b-5dc0-b6b1-59888e4c5d15")].version


function get_version(pkg_name::String)
    for (uuid, pkg) in Pkg.dependencies()
        if pkg.name == pkg_name
            return pkg.version
        end
    end
    return "Package not found in current environment"
end

println(get_version("JLD2"))
Pkg.add("JLD2")

# HPC version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="0.5.15"))
using JLD2
# NBPMscape version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="1.11.2"))