### Bootstrap for median TD by region and correlation between median TDs by region

using JLD2
using NBPMscape
using GLM, Statistics
using DataFrames
#using Pkg
#Pkg.add("Bootstrap")
using Bootstrap
using CSV


# Load median TD simulation replicate data
# Load file containing data for 60,000 replicates (the max so far) and psampled = x%
# ICU sampling proportion
psampled = 1.00 # 0.25
# Data is of format sim replicate (rows) x region (cols) x repeat analysis (slices) and values are median TDs
td_by_simrep_by_region_variability = load("scripts/median_TD_by_region/sim_regionentry/variability_check/td_by_simrep_by_region_variability_60000_$(Int(psampled*100)).jld2", "td_by_simrep_by_region_variability")

# Add 'repeat analyses' together
td_by_simrep_by_region_variability_merged = td_by_simrep_by_region_variability[:,:,1]
for s in 2:size(td_by_simrep_by_region_variability,3)
    td_by_simrep_by_region_variability_merged = vcat(td_by_simrep_by_region_variability_merged
                                                     , td_by_simrep_by_region_variability[:,:,s])
end
# Check completed
size(td_by_simrep_by_region_variability_merged)
td_by_simrep_by_region_variability_merged[60001:120000,:] == td_by_simrep_by_region_variability[:,:,2]

# Number of simulation replicates, regions and 'repeat analyses'
n_sim_reps = size(td_by_simrep_by_region_variability,1)
n_regions = size(td_by_simrep_by_region_variability,2)
n_repeats = size(td_by_simrep_by_region_variability,3)

# Number of bootstrap r
n_bootstrap = 1000

# Create vector to hold bootstrap values
bs_resampled_means = Vector{Union{Missing, Float64}}(missing, n_bootstrap)
bs_resampled_medians = Vector{Union{Missing, Float64}}(missing, n_bootstrap)
bs_resampled_means_2 = Vector{Union{Missing, Float64}}(missing, n_bootstrap)
bs_resampled_medians_2 = Vector{Union{Missing, Float64}}(missing, n_bootstrap)

# Create df to hold results
bs_results_df = DataFrame( zeros(8,39), REGKEY.code)

# Loop through regions
for r in 1:n_regions 

    # Fill bootstrap vectors
    for i in 1:n_bootstrap
        # For simulation replicates
        bs_resample = rand( td_by_simrep_by_region_variability[:,r,1], n_sim_reps)
        bs_resampled_means[i] = mean( filter( isfinite, bs_resample ) )
        bs_resampled_medians[i] = median( filter( isfinite, bs_resample ) )
        
        # For simulation replicates x n repeat analyses
        bs_resample_2 = rand( td_by_simrep_by_region_variability_merged[:,r]
                            , length(td_by_simrep_by_region_variability_merged[:,r]))
        bs_resampled_means_2[i] = mean( filter( isfinite, bs_resample_2) )
        bs_resampled_medians_2[i] = median( filter( isfinite, bs_resample_2 ) )
        
    end
    # For simulation replicates
    # mean
    bs_results_df[1,r] = mean(bs_resampled_means)
    bs_results_df[2,r] = var(bs_resampled_means)
    # median     
    bs_results_df[3,r] = median(bs_resampled_medians)
    bs_results_df[4,r] = var(bs_resampled_medians)

    # # For simulation replicates x n repeat analyses
    # mean
    bs_results_df[5,r] = mean(bs_resampled_means_2)
    bs_results_df[6,r] = var(bs_resampled_means_2)
    # median     
    bs_results_df[7,r] = median(bs_resampled_medians_2)
    bs_results_df[8,r] = var(bs_resampled_medians_2)

end

println(bs_results_df)

# Add column and row names (and transpose)
### Dataframe with regions inc names as rows and columns: median, var, mean, var
bootstrap_results_median_TD = permutedims(bs_results_df)
rename!(bootstrap_results_median_TD, [ "$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_median_TD"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_median_TD_var"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_median_TD"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_median_TD_var"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD_var"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_sim_rep_x_$(n_repeats)_repeats_median_TD_var"
])

# Add ITL2 region codes/names as a new column
ITL2_dict = Dict(zip(REGKEY.code, REGKEY.name))
ITL2_names = [ITL2_dict[k] for k in REGKEY.code]
# Merge pairs into strings (e.g., "key:value" = "code:names")
row_names = [string(k, ":", v) for (k, v) in zip(REGKEY.code, ITL2_names)]
# Insert as column
insertcols!(bootstrap_results_median_TD, 1, :ITL2_region => row_names )

# Save to csv file
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_median_TD_$(n_sim_reps)simreps_x_$(n_repeats)_repeats_p$(Int(psampled*100)).csv", bootstrap_results_median_TD)



### Repeat but using bootstrap for rows, so sampling all regions at the same 
### time for individual simulation replicates. This is because the times to detection (TD) 
### for each region are not independent of each other.

# Create arrays to hold bootstrap values
bs_resampled_means = Array{Union{Missing, Float64}}(missing, n_bootstrap, n_regions)
bs_resampled_medians = Array{Union{Missing, Float64}}(missing, n_bootstrap, n_regions)
bs_resampled_means_2 = Array{Union{Missing, Float64}}(missing, n_bootstrap, n_regions)
bs_resampled_medians_2 = Array{Union{Missing, Float64}}(missing, n_bootstrap, n_regions)

# Create df to hold results
bs_results_df = DataFrame( zeros(8,39), REGKEY.code)

# Fill bootstrap arrays
for i in 1:n_bootstrap

    # Bootstrap sample of simulation replicates 
    # (note that this samples the whole row, i.e. the median TD for each region because
    # they are not independent of each other)
    bs_simrep_row_indices = rand( 1:size(td_by_simrep_by_region_variability,1), n_sim_reps)
    bs_simreps = td_by_simrep_by_region_variability[bs_simrep_row_indices,:,1]

    # Bootstrap sim reps x 'repeat analyses'
    bs_simrep_row_indices_2 = rand( 1:size(td_by_simrep_by_region_variability_merged,1), n_sim_reps*n_repeats)
    bs_simreps_2 = td_by_simrep_by_region_variability_merged[ bs_simrep_row_indices_2 , : , 1 ]

    # Loop through regions
    for r in 1:n_regions 

        # For simulation replicates
        bs_resampled_means[i,r] = mean( filter( isfinite, bs_simreps[:,r] ) )
        bs_resampled_medians[i,r] = median( filter( isfinite, bs_simreps[:,r] ) )
        
        # For simulation replicates x n repeat analyses
        bs_resampled_means_2[i,r] = mean( filter( isfinite, bs_simreps_2[:,r]) )
        bs_resampled_medians_2[i,r] = median( filter( isfinite, bs_simreps_2[:,r] ) )
        
    end

end

# Loop through regions again to gather results from bootstrap samples
for r in 1:n_regions
    
    # Fill in resluts df
    # For simulation replicates
    # mean
    bs_results_df[1,r] = mean(bs_resampled_means[:,r])
    bs_results_df[2,r] = var(bs_resampled_means[:,r])
    # median     
    bs_results_df[3,r] = median(bs_resampled_medians[:,r])
    bs_results_df[4,r] = var(bs_resampled_medians[:,r])

    # # For simulation replicates x n repeat analyses
    # mean
    bs_results_df[5,r] = mean(bs_resampled_means_2[:,r])
    bs_results_df[6,r] = var(bs_resampled_means_2[:,r])
    # median     
    bs_results_df[7,r] = median(bs_resampled_medians_2[:,r])
    bs_results_df[8,r] = var(bs_resampled_medians_2[:,r])
end

println(bs_results_df)

# Add column and row names (and transpose)
### Dataframe with regions inc names as rows and columns: median, var, mean, var
bootstrap_by_simrep_results_median_TD = permutedims(bs_results_df)
rename!(bootstrap_by_simrep_results_median_TD, [ "$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_median_TD"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_median_TD_var"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_median_TD"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_median_TD_var"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD"
                            ,"$(n_bootstrap)_bs_mean_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD_var"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_simreps_x_$(n_repeats)_repeats_median_TD"
                            ,"$(n_bootstrap)_bs_median_$(n_sim_reps)_sim_rep_x_$(n_repeats)_repeats_median_TD_var"
])

# Add ITL2 region codes/names as a new column
ITL2_dict = Dict(zip(REGKEY.code, REGKEY.name))
ITL2_names = [ITL2_dict[k] for k in REGKEY.code]
# Merge pairs into strings (e.g., "key:value" = "code:names")
row_names = [string(k, ":", v) for (k, v) in zip(REGKEY.code, ITL2_names)]
# Insert as column
insertcols!(bootstrap_by_simrep_results_median_TD, 1, :ITL2_region => row_names )

# Save to csv file
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_by_sim_rep_results_median_TD_$(n_sim_reps)simreps_x_$(n_repeats)_repeats_p$(Int(psampled*100)).csv", bootstrap_by_simrep_results_median_TD)


### Repeat with Bootstrap package
bs_mean = bootstrap(mean, filter(isfinite,td_by_simrep_by_region_variability[:,1,1]), BasicSampling(n_bootstrap))
bs_median = bootstrap(median, filter(isfinite,td_by_simrep_by_region_variability[:,1,1]), BasicSampling(n_bootstrap))


### Bootstrap correlation between regions

### Because the correlation values are not independent, the bootstrapping needs to 
# be on a simulation replicate level 
# TODO - The same is probably true for the median TD as well above

# Number of bootstraps
n_bootstrap = 1000

### Computing correlation between ITL2 regions for TDs (for all simulation replicates)
n = size(td_by_simrep_by_region_variability, 2)
corrmat     = Matrix{Union{Missing, Float64}}(undef, n, n)
corrmat_idx = Matrix{Union{Missing, Float64}}(undef, n, n)
corrmat_bs_array = Array{Union{Missing, Float64}}(undef, n, n, n_bootstrap)
corrmat2     = Matrix{Union{Missing, Float64}}(undef, n, n)
corrmat_simreps_x_repeats_bs_array = Array{Union{Missing, Float64}}(undef, n, n, n_bootstrap)

for bs in 1:n_bootstrap
    # Bootstrap sample of simulation replicates 
    # (note that this samples the whole row, i.e. the median TD for each region because
    # they are not independent of each other)
    bs_simrep_row_indices = rand( 1:size(td_by_simrep_by_region_variability,1), n_sim_reps)
    bs_simreps = td_by_simrep_by_region_variability[bs_simrep_row_indices,:,1]

    # Bootstrap sample of simulation replicates and repeat analyses
    bs_simrep_x_repeats_row_indices = rand( 1:size(td_by_simrep_by_region_variability_merged,1), n_sim_reps*n_repeats)
    bs_simreps_x_repeats = td_by_simrep_by_region_variability_merged[bs_simrep_x_repeats_row_indices,:,1]

    for i in 1:n
        for j in i:n
            ### For simulation replicates
            # Extract columns as vectors
            x = bs_simreps[:, i, 1]
            y = bs_simreps[:, j, 1]
            # Find indices where both are not missing...
            #idx = .!ismissing.(x) .& .!ismissing.(y)
            # ... and are not infinite (i.e. there is no ICU case in that simulation replicate)
            idx = isfinite.(x) .& isfinite.(y)
            #TEST
            #sum(isfinite.(x))
            #sum(isfinite.(y))
            #sum( isfinite.(x) .& isfinite.(y) )
            # Compute correlation on complete cases
            #corr = cor(skipmissing(x[idx]), skipmissing(y[idx]))
            x_filtered = skipmissing(x[idx])[:]
            y_filtered = skipmissing(y[idx])[:]
            corr = sum(idx) <= 1 ? missing : cor(x_filtered, y_filtered)
            corrmat[i, j] = corr
            corrmat[j, i] = corr  # Symmetric
            #corrmat_idx[i, j] = sum(idx)
            #corrmat_idx[j, i] = sum(idx) # Symmetric

            ### For simulation replciates x repeat analyses
            # Extract columns as vectors
            x2 = bs_simreps_x_repeats[:, i, 1]
            y2 = bs_simreps_x_repeats[:, j, 1]
            # Find indices where both are not missing...
            #idx = .!ismissing.(x) .& .!ismissing.(y)
            # ... and are not infinite (i.e. there is no ICU case in that simulation replicate)
            idx2 = isfinite.(x2) .& isfinite.(y2)
            # Compute correlation on complete cases
            x_filtered2 = skipmissing(x2[idx2])[:]
            y_filtered2 = skipmissing(y2[idx2])[:]
            corr2 = sum(idx2) <= 1 ? missing : cor(x_filtered2, y_filtered2)
            corrmat2[i, j] = corr2
            corrmat2[j, i] = corr2  # Symmetric
            
        end
    end
    corrmat_bs_array[:,:,bs] = corrmat
    corrmat_simreps_x_repeats_bs_array[:,:,bs] = corrmat2
end
#println(corrmat_idx)
#mean(corrmat_idx)

bs_corrmat_mean   = mean(corrmat_bs_array, dims=3)
bs_corrmat_median = median(corrmat_bs_array, dims=3)
bs_corrmat_var    = var(corrmat_bs_array, dims=3)

bs_corrmat_mean_df = DataFrame( dropdims( bs_corrmat_mean, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_mean_df, 1, :ITL2_region => row_names )

bs_corrmat_median_df = DataFrame( dropdims(bs_corrmat_median, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_median_df, 1, :ITL2_region => row_names )

bs_corrmat_var_df = DataFrame( dropdims(bs_corrmat_var, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_var_df, 1, :ITL2_region => row_names )

# Save to csv file
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_mean_corr_median_TD_$(n_sim_reps)simreps_p$(Int(psampled*100)).csv"
            , bs_corrmat_mean_df)
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_median_corr_median_TD_$(n_sim_reps)simreps_p$(Int(psampled*100)).csv"
            , bs_corrmat_median_df)
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_var_corr_median_TD_$(n_sim_reps)simreps_p$(Int(psampled*100)).csv"
            , bs_corrmat_var_df)

bs_corrmat_simreps_x_repeats_mean   = mean(corrmat_simreps_x_repeats_bs_array, dims=3)
bs_corrmat_simreps_x_repeats_median = median(corrmat_simreps_x_repeats_bs_array, dims=3)
bs_corrmat_simreps_x_repeats_var    = var(corrmat_simreps_x_repeats_bs_array, dims=3)

bs_corrmat_simreps_x_repeats_mean_df = DataFrame( dropdims(bs_corrmat_simreps_x_repeats_mean, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_simreps_x_repeats_mean_df, 1, :ITL2_region => row_names )

bs_corrmat_simreps_x_repeats_median_df = DataFrame( dropdims(bs_corrmat_simreps_x_repeats_median, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_simreps_x_repeats_median_df, 1, :ITL2_region => row_names )

bs_corrmat_simreps_x_repeats_var_df = DataFrame( dropdims(bs_corrmat_simreps_x_repeats_var, dims=3), Symbol.(REGKEY.code)) #names(td_by_simrep_by_region)))
insertcols!(bs_corrmat_simreps_x_repeats_var_df, 1, :ITL2_region => row_names )

# Save to csv file
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_mean_corr_median_TD_$(n_sim_reps)simreps_x_$(n_repeats)_p$(Int(psampled*100)).csv"
            , bs_corrmat_simreps_x_repeats_mean_df)
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_median_corr_median_TD_$(n_sim_reps)simreps_x_$(n_repeats)_p$(Int(psampled*100)).csv"
            , bs_corrmat_simreps_x_repeats_median_df)
CSV.write("scripts/median_TD_by_region/sim_regionentry/bootstrap$(n_bootstrap)_results_var_corr_median_TD_$(n_sim_reps)simreps_x_$(n_repeats)_p$(Int(psampled*100)).csv"
            , bs_corrmat_simreps_x_repeats_var_df)