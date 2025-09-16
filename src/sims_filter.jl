"""
Function: sims_filter

Filters file containing simulations output by {simtree} or {simforest} for the G dataframe
and then rows containing values for ticu (infected individuals admitted to ICU) and 
tgp (infected individuals that presented to GP).
The main reason for doing this filtering is to reduce file size, which is useful for storage,
but more importantly, very large .jld2 files are difficult to reload.

# Arguments
'sims_file':        File containing simulations output by {simtree} or {simforest}
'sims_object_name': The name of the object saved to the .jld2 file. Without this name the data cannot 
                    be loaded from the .jld2 file
'filter_variables': Currently this must be tgp and ticu formatted as a 2-element vector ["tgp","ticu"]

# Returns
A dataframe in the same format as the G dataframe but with the rows filtered for those with non-infinite values in the tgp and ticu columns:
 
 Row │ pid                                tinf     tgp      thospital  ticu      trecovered  severity    iscommuter  homeregion  commuteregion  generation  F      G      H           infector_age  infectee_age  simid                             
     │ String                             Float64  Float64  Float64    Float64   Float64     Symbol      Bool        String      String         Int64       Int64  Int64  Float64     Int8?         Int8          String
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ f6f6b432-8e44-11f0-3306-27246943…  19.0561  23.6572    25.5282  Inf       31.7309     severe      true        TLE4        TLE4           2           3      0      11.8131      84            90           f6f6b432-8e44-11f0-3306-27246943…
   2 │ f6f6b432-8e44-11f0-3306-27246943…  18.5838  20.7196   Inf       Inf       Inf         moderate    true        TLE4        TLE4           3           8      39     0.411181     38            49           f6f6b432-8e44-11f0-3306-27246943…
  ...| ...                                ...      ...      ...        ...       ...         ...         ...         ...         ...            ...         ...

# Example
    # Run filter function
    [sims_G_icu_filter,sims_G_gp_filter] = sims_filter( sims_file = "sims.jld2", sims_object_name = "sims", filter_variables = ["tgp","ticu"])
    # Save output filtered data to file
    @save "sims_filtered_G_icu.jld2" sims_G_icu_filter
    @save "sims_filtered_G_gp.jld2" sims_G_gp_filter

"""
using DataFrames
using JLD2

function sims_filter(; sims_file = "sims.jld2"
                    , sims_object_name = "sims"
                    , filter_variables = ["tgp","ticu"])

    # Load file
    #sims = load("covidlike-1.1.1-sims.jld2", "sims")
    sims = load(sims_file, sims_object_name)

    # Create vector to store filtered dataframes
    sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]
    sims_G_gp_filter  = [DataFrame() for _ in 1:length(sims)]

    # Loop through replicates
    for s in 1:length(sims)
        
        try
            fo = sims[s]
            # filter for G, which is the dataframe containing infection information,
            # and only retain cases (rows) that:
            # - progressed to ICU (i.e. capable of detection under ICU sampling methodology)
            G_icu = fo.G[ isfinite.(fo.G.ticu), : ]
            # - progressed to GP (i.e. capable of detection under primarycare sampling methodology)
            G_gp = fo.G[ isfinite.(fo.G.tgp), : ]
            
            # If the filtered df has rows then add them to the vector holding the dfs,
            # otherwise add 'missing' to the vector
            size(G_icu,1) > 0 ? sims_G_icu_filter[s] = G_icu : next #sims_G_icu_filter[s] = missing
            size(G_gp, 1) > 0 ? sims_G_gp_filter[s] = G_gp  : next #sims_G_gp_filter[s]  = missing

        catch err
            @warn "Error in iteration" exception = err
            println("Error in sim number: ",s)
            continue
        end     

    end
    
    return( [sims_G_icu_filter, sims_G_gp_filter] )

end