# Read in multiple simulation files and combine for analysis of a larger group of simulation replicates

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

sim_root_path = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/nrep_1000/sims/"

sim_files_1 = readdir("$(sim_root_path)779103/"; join=true)
sim_files_2 = readdir("$(sim_root_path)779443/"; join=true)
sim_files_3 = "covidlike-1.0-sims-regionentry_filtered_G_icu.jld2"
sim_files_4 = "covidlike-1.0-sims-regionentry-nrep5000_filtered_G_icu.jld2"
sim_files_5 = "covidlike-1.0-sims-regionentry-nrep5000_2_filtered_G_icu.jld2"
# 2nd batch of files to increase size of sim reps
sim_files_6 = readdir("$(sim_root_path)787041/"; join=true)
sim_files_7 = readdir("$(sim_root_path)787313/"; join=true)
sim_files_8 = readdir("$(sim_root_path)787315/"; join=true)
sim_files_9 = readdir("$(sim_root_path)787353/"; join=true)

sim_files = vcat(sim_files_1, sim_files_2, sim_files_3, sim_files_4, sim_files_5)

#nrep = 1000
#sims = Vector{Any}(undef, nrep * length(sim_files))

for i in 1:length(sim_files)
    println("Processing ", sim_files[i])
    # Load simulation replicates (G filtereed for ICU cases only) from file
    sims_temp = load(sim_files[i], "sims_G_icu_filter")

    # Number of replicates 
    nrep = length(sims_temp)
    
    # Append sim replicates from each file to a single vector with a df for each sim replicate
    if i == 1
        sims_G_icu_filter = sims_temp
    else
        sims_G_icu_filter = append!(sims_G_icu_filter, sims_temp)
    end

end

# Save combined object to file
@save "covidlike-1.0-sims-regionentry_filtered_G_icu_combined.jld2" sims_G_icu_filter


### 2nd batch of files to increase size of sim reps
sim_files_6 = readdir("$(sim_root_path)787041/"; join=true)
sim_files_7 = readdir("$(sim_root_path)787313/"; join=true)
sim_files_8 = readdir("$(sim_root_path)787315/"; join=true)
sim_files_9 = readdir("$(sim_root_path)787353/"; join=true)

sim_files_2 = vcat(sim_files_6, sim_files_7, sim_files_8, sim_files_9)

# Check current length of file in terms of sim reps
length(sims_G_icu_filter)

for j in 1:length(sim_files_2)
    println("Processing ", sim_files_2[j])
    # Load simulation replicates (G filtereed for ICU cases only) from file
    sims_temp = load(sim_files_2[j], "sims_G_icu_filter")

    # Number of replicates 
    nrep = length(sims_temp)
    
    # Append sim replicates from each file to a single vector with a df for each sim replicate
    #if i == 1
    #    sims_G_icu_filter = sims_temp
    #else
        sims_G_icu_filter = append!(sims_G_icu_filter, sims_temp)
    #end

end

# Check length of file in terms of sim reps now
length(sims_G_icu_filter)

# Save combined object to file
@save "covidlike-1.0-sims-regionentry_filtered_G_icu_combined_nrep60000.jld2" sims_G_icu_filter

# Check loads
sims_G_icu_filter = load("covidlike-1.0-sims-regionentry_filtered_G_icu_combined_nrep60000.jld2", "sims_G_icu_filter")

##########################################
# HPC version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="0.5.15"))
using JLD2
# NBPMscape version of JLD2
using Pkg
Pkg.add(Pkg.PackageSpec(name="JLD2", version="1.11.2"))