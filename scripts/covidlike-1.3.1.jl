#= 
- Simulate multiple replicates using default covid-like parameters 
- Plus age disaggregation of contact numbers and assortativity
- Plus age disaggregation of commuting (only people aged 16 years and above)
- Plus age disaggregation of infection severity
- Plus more detailed care pathways and rates
- Plus death COMPARTMENT
- Plus recording whether infection was imported
- Prep for sub-sampling
- Serialise 
=#

using Pkg
Pkg.instantiate()
Pkg.resolve()

using NBPMscape 
using JLD2 

NREPS = 1000 #1000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100), 1:NREPS )
@save "covidlike-1.3.1-sims-nrep1000.jld2" sims

#TODO NOT WORKING YET
for i in 1:length(sims)
    sim_temp = sims[i]
    #var_name = Symbol("sims_$(i)")   # create a symbol like :df_1, :df_2, ...
    #eval(:( $(var_name) = $sim_temp ))  
    @save "covidlike-1.1.1-sims/covidlike-1.1.1-sims_$i.jld2" sim_temp #$(var_name)
end

# Filter for G and ICU and GP cases and save (reduces file size and more likely to be able to reload)
# Create vector to store filtered dataframes
sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]
sims_G_gp_filter  = [DataFrame() for _ in 1:length(sims)]

# Loop through replicates
for s in 1:length(sims)
    
    #println("file: ",sfn,", sim number: ",s)
    try
        fo = sims[s]
        # filter for G, which is the dataframe containing infection information,
        # and only retain cases (rows) that progressed to ICU (i.e. capable of detection under ICU sampling methodology)
        G_icu = fo.G[ isfinite.(fo.G.ticu), : ]
        G_gp = fo.G[ isfinite.(fo.G.tgp), : ]

        if size(G_icu,1) > 0
            sims_G_icu_filter[s] = G_icu
        else
            sims_G_icu_filter[s] = missing
            #continue        
        end

        if size(G_gp,1) > 0
            sims_G_gp_filter[s] = G_gp
        else
            sims_G_gp_filter[s] = missing
            #continue        
        end

    catch err
        @warn "Error in iteration" exception = err
        println("Error in sim number: ",s)
        continue
    end     

end
# Save filtered data to file
@save "covidlike-1.3.1-sims_filtered_G_icu_nrep10.jld2" sims_G_icu_filter
@save "covidlike-1.3.1-sims_filtered_G_gp_nrep10.jld2" sims_G_gp_filter

# Also save other constituent parts of sims