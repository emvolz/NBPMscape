# Read in multiple simulation files and combine for analysis of a larger group of
# simulation replicates

using JLD2
using DataFrames
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

#sim_root_path = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/"
sim_root_path_1 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/"
sim_root_path_2 = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/951401/"

sim_files_gp_1 = readdir("$(sim_root_path_1)G_filtered_gp/"; join=true)
sim_files_icu_1 = readdir("$(sim_root_path_1)G_filtered_icu/"; join=true)

#sim_files_gp_2 = readdir("$(sim_root_path_2)G_filtered_gp/"; join=true)
#sim_files_icu_2 = readdir("$(sim_root_path_2)G_filtered_icu/"; join=true)

sim_files_icu = vcat(sim_files_icu_1, sim_files_icu_2)
sim_files_gp = vcat(sim_files_gp_1, sim_files_gp_2)

# If only a single folder containing files
sim_files_icu = sim_files_icu_1
sim_files_gp = sim_files_gp_1

nrep = 1000
#sims = Vector{Any}(undef, nrep * length(sim_files))

# Create vector to store filtered dataframes
sims_G_icu_filter = [DataFrame() for _ in 1:length(nrep)]
sims_G_gp_filter  = [DataFrame() for _ in 1:length(nrep)]

# Add all simulation replicates from multiple files to a single object 
for i in 1:length(sim_files_icu)
    println("Processing ", sim_files_icu[i])
    # Load simulation replicates (G filtereed for ICU cases only) from file
    sims_temp = load(sim_files_icu[i], "sims_G_icu_filter")

    # Number of replicates 
    nrep = length(sims_temp)
    
    # Append sim replicates from each file to a single vector with a df for each sim replicate
    if i == 1
        sims_G_icu_filter = sims_temp
    else
        sims_G_icu_filter = append!(sims_G_icu_filter, sims_temp)
    end

end

# Check length of file in terms of sim reps now
nrep_sims_G_icu_filter = length(sims_G_icu_filter)

# Save combined object to file
#@save "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2" sims_G_icu_filter
#@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2" sims_G_icu_filter
@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter)_955898.jld2" sims_G_icu_filter

# Add all simulation replicates from multiple files to a single object
for i in 1:length(sim_files_gp)
    println("Processing ", sim_files_gp[i])
    # Load simulation replicates (G filtereed for GP cases only) from file
    sims_temp = load(sim_files_gp[i], "sims_G_gp_filter")

    # Number of replicates 
    nrep = length(sims_temp)
    
    # Append sim replicates from each file to a single vector with a df for each sim replicate
    if i == 1
        sims_G_gp_filter = sims_temp
    else
        sims_G_gp_filter = append!(sims_G_gp_filter, sims_temp)
    end

end

# Check length of file in terms of sim reps now
nrep_sims_G_gp_filter = length(sims_G_gp_filter)

# Save combined object to file
#@save "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2" sims_G_gp_filter
#@save "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2" sims_G_gp_filter
@save "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter)_955898.jld2" sims_G_gp_filter

# Check loads
#sims_G_icu_filter = load("covidlike-1.1.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2", "sims_G_icu_filter")
#sims_G_gp_filter = load("covidlike-1.1.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2", "sims_G_gp_filter")
sims_G_icu_filter = load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2", "sims_G_icu_filter")
sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2", "sims_G_gp_filter")

# Count how many simulations contain results
function n_sims_w_results(; sims_file_icu, sims_file_gp)
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

n_sims_w_results( sims_file_icu = "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep53000.jld2" 
                , sims_file_gp = "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep53000.jld2") 
# ICU: 90% and GP: 90% 

n_sims_w_results( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep84000.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep84000.jld2") 
# Combinatin of HPC run 851401+854085 
# ICU: 35% and GP: 34% (because R was too low at ~1)

n_sims_w_results( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2") 
# ICU: 80% and GP: 80%


##########################################
# HPC version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="0.5.15"))
using JLD2
# NBPMscape version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="1.11.2"))