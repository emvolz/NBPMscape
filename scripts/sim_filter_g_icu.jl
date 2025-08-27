# Script to reduce file size of simulation

using JLD2
using DataFrames

# load simulation
sim_file_name_root_1 = "covidlike-1.0-sims-regionentry"
sim_file_name_root_2 = "covidlike-1.0-sims-regionentry-nrep5000" # sim replicate 4396 appears to be empty
sim_file_name_root_3 = "covidlike-1.0-sims-regionentry-nrep5000_2" # sim replicate 2333 appears to be empty

sim_file_names = [ sim_file_name_root_1, sim_file_name_root_2, sim_file_name_root_3]

for sfn in sim_file_names
    sims = load("$(sfn).jld2" , "sims")

    # Create vector to store filtered dataframes
    sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]

    # Loop through replicates
    for s in 1:length(sims)
    
        #println("file: ",sfn,", sim number: ",s)
        try
            fo = sims[s]
            # filter for G, which is the dataframe containing infection information,
            # and only retain cases (rows) that progressed to ICU (i.e. capable of detection under ICU sampling methodology)
            G = fo.G[ isfinite.(fo.G.ticu), : ]

            if size(G,1) > 0
                sims_G_icu_filter[s] = G
            else
                continue        
            end

        catch err
            @warn "Error in iteration" exception = err
            println("Error in file: ",sfn,", sim number: ",s)
            continue
        end     

    end
# Save filtered data to file
@save "$(sfn)_filtered_G_icu.jld2" sims_G_icu_filter

sims_G_icu_filter = missing

end


# Check can reload
sims_G_icu_filter = load("covidlike-1.0-sims-regionentry_filtered_G_icu.jld2", "sims_G_icu_filter")
obj_type = typeof(sims_G_icu_filter)
obj_type.parameters
obj_type.name.wrapper
Base.namedtupletypeparameters(obj_type)