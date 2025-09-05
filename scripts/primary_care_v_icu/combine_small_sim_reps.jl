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

sim_root_path = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/"
#sim_root_path = "~/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/"

sim_files_gp = readdir("$(sim_root_path)G_filtered_gp/"; join=true)
sim_files_icu = readdir("$(sim_root_path)G_filtered_icu/"; join=true)

#sim_files = vcat(sim_files_1, sim_files_2)

#nrep = 1000
#sims = Vector{Any}(undef, nrep * length(sim_files))

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
@save "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2" sims_G_icu_filter

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
@save "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2" sims_G_gp_filter

# Check loads
sims_G_icu_filter = load("covidlike-1.1.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2", "sims_G_icu_filter")
sims_G_gp_filter = load("covidlike-1.1.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2", "sims_G_gp_filter")

##########################################
# HPC version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="0.5.15"))
using JLD2
# NBPMscape version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="1.11.2"))