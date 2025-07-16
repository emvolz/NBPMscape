### Compute median time to detection (TD) for all ITL2 regions individually

#=
Outline
 1) Extract event timing information by region from simulation 
 2) Fit statistical distributions to event timing information 
 3) Report estimated TD by region including plots and df
 
Description
 Extract region specific information from each forest in the outbreak simulation. 
 However, not all forests will contain ICU infections (able to be sampled) in a 
 particular ITL2 region. In that case, a robust way of estimating it would be to fit a
 parametric distribution (e.g. gamma) to the TD's. 
 The simulation has maximum run time and so truncation needs to be accounted for.
 Accounting for truncation is more difficult for treport and so better to use tinf as 
 this is the parameter that is actually truncated by maxtime.
 Possible to fit to tinf and then fit separately to treport - tinf, which is not truncated.
 These can then be summed to obtain an estimate of treport (= tinf + (treport - tinf))
 Finally report median estimates for TD (earliest treport) by region.
=#

using JLD2
using NBPMscape
using GLM, Statistics
using Distributions
using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using DataFrames
import UUIDs 
import StatsBase
using StatsBase
using Interpolations
import SpecialFunctions as SF 
using Plots 
using LinearAlgebra
using Optim
using RData 
using CSV 
using PowerAnalyses

# Load simulation
#sims = load("covidlike-1.0-sims.jld2" , "sims") # preliminary simulation
sims = load("covidlike-1.0-sims-regionentry.jld2" , "sims") # preliminary simulation re-run with regionentry adjusted

# ITL2 region list
REGKEY.code

#Testing
#for s in 1:length(sims)
#    G = sims[s][1]
#    println(names(G))
#    D = sims[s][2]
#    n_imports = sims[s][3]
#    t_imports = sims[s][4]
#println(G)
#end

#fo = sims[1]

# Sampling proportion within region
psampled = 0.1 #1.0 #0.1
# Time to process sample and report results / declare detection
turnaroundtime = 3

#length(REGKEY.code)
#s=1

# Create dictionary to hold the simulated time of infection for each ITL2 region
tinf_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
tsample_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
tinf_tsample_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
tinf_treport_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
treport_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
# Create df to hold the number of ICU cases by region in each simulation replicate
n_icu_by_region_by_simrep = DataFrame( zeros( length(sims) , length(REGKEY.code)), REGKEY.code )
# Create df to record the TD by simulation replicate and by region
td_by_simrep_by_region = DataFrame( [fill(Inf, length(sims)) for _ in 1:length(REGKEY.code)], REGKEY.code)

# Loop through each simulation replicate and compile information
for s in 1:length(sims)
    #s=1
    fo = sims[s]
    
    # Filter for ICU cases 
    G = fo.G[ isfinite.(fo.G.ticu), : ]
    #println( names(G) )
    
    # Record number of ICU cases by region and by simulation
    icu_count_by_region = countmap(G.homeregion)
    
    for (k, v) in icu_count_by_region
        n_icu_by_region_by_simrep[ s, k ] = v
    end
    #println(icu_count_by_region)
    #println(n_icu_by_region_by_simrep[1,:])

    # List regions with an ICU admission (i.e. possibility of getting a sample and detection)
    icu_regions = unique(G.homeregion)
    
    # Loop through each region and record time of infection stats
    for r in icu_regions #REGKEY.code
        #Test
        #r=icu_regions[7]
        
        g_region = G[ (G.homeregion).==r, :]
        
        # Sample size 
	    n = rand( Binomial( size(g_region,1), psampled) ) #p.psampled ) )
	    # Subsample of ICU cases within region
	    g_region_sub = g_region[sample( 1:size(g_region,1), n, replace=false ), :]

        # Check if there are any rows before continuing
        if size(g_region_sub,1) == 0 
            continue
        elseif size(g_region_sub,1) > 0 
            
            #= Sample time has uniform distribution between time of admission to ICU and time of recovery
               TODO THIS MAY NEED TO BE UPDATED FOR LATER VERSIONS WITH MORE COMPLEX CARE PATHWAYS
               NOTE TIME BETWEEN ticu and trecovered can be large and so with uniform distribution of 
               sampling the tsample can be a long time after tinf. 
               Example I saw was tinf = 48.9, ticu = 52.6, trecovered = 91.1, and treport = 90.7
               If we're modelling 100% of ICU cases being sampled then more likely to sampled closer to ticu
            =#
            
            # Generate sample times
            #tsample = (size(g_region,1)>0) ? map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region) ) : []
            tsample = map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region_sub) )
            g_region_sub.tsample = tsample 
            #println(g_region_sub[:,[2,3,4,5,6,16]])

            # Simulate reports times
            #treport = (size(g_region,1)>0) ? (tsample .+ turnaroundtime) : []
            treport = (tsample .+ turnaroundtime)
            g_region_sub.treport = treport
            #println(g_region_sub[:,[2,3,4,5,6,16,17]])

            # Find infection (row) with minimum report (detection) time
            min_treport, min_treport_row_index = findmin( g_region_sub.treport )

            ## Add times to respective dictionaries
            # Times of infection for ICU admissions
            push!(tinf_by_region_dict[ r ], g_region_sub.tinf[min_treport_row_index] )
            # Times of sampling for ICU admissions
            push!(tsample_by_region_dict[ r ], g_region_sub.tsample[min_treport_row_index] )
            # Times of report / detection for ICU admissions
            push!(treport_by_region_dict[ r ], g_region_sub.treport[min_treport_row_index] )
            # Time between infection and sampling for ICU admissions
            push!(tinf_tsample_by_region_dict[ r ], g_region_sub.tsample[min_treport_row_index] - g_region_sub.tinf[min_treport_row_index] )
            # Time between infection and reporting of sample results for ICU admissions
            push!(tinf_treport_by_region_dict[ r ], g_region_sub.treport[min_treport_row_index] - g_region_sub.tinf[min_treport_row_index] )
            
            # Record TD (= min(treport)) by simulation replicate and by region
            td_by_simrep_by_region[ s , r ] = min_treport
        end
                
    end

end

# Check how many tinf = 0 for each region (this is the tinf relating to the earliest treport)
tinf_zero_count_by_region = DataFrame( zeros( 1, length( tinf_by_region_dict ) ), REGKEY.code )
tinf_zero_perc_by_region = DataFrame( zeros( 1, length( tinf_by_region_dict ) ), REGKEY.code )
for (k, v) in tinf_by_region_dict
    tinf_zero_count_by_region[ 1, k ] = count(==(0), v)
    tinf_zero_perc_by_region[ 1, k ] = 100 *count(==(0), v) / length(v)
end
println(tinf_zero_count_by_region)
minimum(tinf_zero_count_by_region[1,:])
maximum(tinf_zero_count_by_region[1,:])
println(tinf_zero_perc_by_region)
minimum(tinf_zero_perc_by_region[1,:])
maximum(tinf_zero_perc_by_region[1,:])


# Check how many treport = 0 for each region
treport_zero_count_by_region = DataFrame( zeros( 1, length( treport_by_region_dict ) ), REGKEY.code )
treport_zero_perc_by_region = DataFrame( zeros( 1, length( treport_by_region_dict ) ), REGKEY.code )
for (k, v) in treport_by_region_dict
    treport_zero_count_by_region[ 1, k ] = count(==(0), v)
    treport_zero_perc_by_region[ 1, k ] = 100 *count(==(0), v) / length(v)
end
println( treport_zero_count_by_region )
minimum(treport_zero_count_by_region[1,:])
maximum(treport_zero_count_by_region[1,:])
println( treport_zero_perc_by_region )
minimum(treport_zero_perc_by_region[1,:])
maximum(treport_zero_perc_by_region[1,:])

# Check how many treport - tinf = 0 for each region
tinf_treport_zero_count_by_region = DataFrame( zeros( 1, length( tinf_treport_by_region_dict ) ), REGKEY.code )
tinf_treport_zero_perc_by_region = DataFrame( zeros( 1, length( tinf_treport_by_region_dict ) ), REGKEY.code )
for (k, v) in tinf_treport_by_region_dict
    tinf_treport_zero_count_by_region[ 1, k ] = count(==(0), v)
    tinf_treport_zero_perc_by_region[ 1, k ] = 100 * count(==(0), v) / length(v)
end
println( tinf_treport_zero_count_by_region )
minimum( tinf_treport_zero_count_by_region[1,:] )
maximum( tinf_treport_zero_count_by_region[1,:] )
println( tinf_treport_zero_perc_by_region )
minimum( tinf_treport_zero_perc_by_region[1,:] )
maximum( tinf_treport_zero_perc_by_region[1,:] )

# Compute % of simulation replicates that include an ICU case by region
println( n_icu_by_region_by_simrep )
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/n_icu_cases_by_region_by_simrep.csv", n_icu_by_region_by_simrep)
perc_simreps_icu_by_region = DataFrame( zeros(1,39), REGKEY.code)
for i in 1:ncol( n_icu_by_region_by_simrep )
    perc_simreps_icu_by_region[1,i] = 100 * sum( n_icu_by_region_by_simrep[:,i] .> 0 ) / nrow(n_icu_by_region_by_simrep)
end
println( perc_simreps_icu_by_region )
minimum(perc_simreps_icu_by_region[1,:])
maximum(perc_simreps_icu_by_region[1,:])

names(perc_simreps_icu_by_region) == names(tinf_zero_perc_by_region)
scatter(Vector(perc_simreps_icu_by_region[1,:])
        , Vector(tinf_zero_perc_by_region[1,:])
        , xlabel="% of sim replicates with ICU cases", ylabel = "% of TD with tinf = 0"
        , title = "By ITL2 region"
        , label = "ITL2 regions")

# Compute median TD by region across all simulation replicates...
td_median_by_region = DataFrame( zeros(1,39), REGKEY.code)
for i in 1:ncol( td_median_by_region )
    td_median_by_region[1,i] = median( filter( isfinite, td_by_simrep_by_region[:,i] ) )
end
# ... and median TD by simulation replicate across all regions
td_median_by_simrep = Vector{Union{Missing, Float64}}(missing, length(sims)) #DataFrame( zeros(1000,1), :auto )
td_min_by_simrep = Vector{Union{Missing, Float64}}(missing, length(sims)) #DataFrame( zeros(1000,1), :auto )
for j in 1:length( td_median_by_simrep )
    #td_median_by_simrep[j,1] = median( filter( isfinite, td_by_simrep_by_region[j,:] ) )
    data_inf_filtered = filter( isfinite, [ getindex( td_by_simrep_by_region, j, k) for k in 1:size( td_by_simrep_by_region, 2)] ) 
    td_median_by_simrep[j] = isempty( data_inf_filtered ) ? missing : median( data_inf_filtered )
    td_min_by_simrep[j] = isempty( data_inf_filtered ) ? missing : minimum( data_inf_filtered )
end
#TEST
#histogram(data_inf_filtered, breaks = 20)

# Plot to compare the median TD values by region and by simulation replicate
median_td_by_region = [getindex( td_median_by_region, 1, i) for i in 1:size( td_median_by_region, 2)] 
histogram( median_td_by_region; #bins=30
              normalize=:pdf
             , label="median TDs by region (all sim reps)"
             , color=:blue
             , alpha=0.5
             , legend=:topleft
             , title="Comparison of median TDs by region and by simulation replicate" #"$(r) $(REGKEY[REGKEY.code .== r, :name][1])" 
             , xlabel="Days", ylabel = "Density"
             , size = (800, 400)
             , xlimit = [0,100]
             , ylimit = 0.2
             , bins = 0:2:100
)
#median_td_by_simrep = [getindex( td_median_by_simrep, 1, i) for i in 1:size( td_median_by_simrep, 1)] 
histogram!( td_median_by_simrep; #bins=50
              normalize=:pdf
             , label="median TDs by sim rep (all regions)"
             , color=:red, alpha=0.5#, legend=:topleft
             , bins = 0:2:100
             )
vline!( [median(median_td_by_region)], color=:blue
            , label="median of median TDs by region (all sim reps)"
            , linewidth = 3) 
vline!( [median(skipmissing(td_median_by_simrep))], color=:red
        , label="median of median TDs by sim rep (all regions)", linewidth = 3  )

# Add minimum TD values by region for each simulation replicate
histogram!( td_min_by_simrep; #bins=50
              normalize=:pdf
             , label="minimum of regional TDs by sim rep"
             , color=:green, alpha=0.5 #, legend=:topleft
             , bins = 0:2:60
             )
vline!( [median(skipmissing(td_min_by_simrep))], color=:green
        , label="median of minimum of regional TDs by sim rep", linewidth = 3  )
#vline!( [quantile( skipmissing(td_min_by_simrep) ,0.025)], color=:green
#        , label="95% CI: minimum of regional TDs by sim rep"
#        , linewidth = 3  , linestyle = :dash)
#vline!( [quantile( skipmissing(td_min_by_simrep) ,0.975)], color=:green
#        , linewidth = 3  , linestyle = :dash)

println("Median of the median TD by region = ", round(median(median_td_by_region), digits=1))
println("Median of the minimum of regional TDs by sim rep = ", round( median(skipmissing(td_min_by_simrep)), digits=1))
println("with 95% CI: ",round(quantile( skipmissing(td_min_by_simrep) ,0.025), digits = 1)
        ," to ", round(quantile( skipmissing(td_min_by_simrep) ,0.975), digits=1))

### Computing correlation between ITL2 regions for TDs (for all simulation replicates)
n = size(td_by_simrep_by_region, 2)
corrmat = Matrix{Union{Missing, Float64}}(undef, n, n)
corrmat_idx = Matrix{Union{Missing, Float64}}(undef, n, n)

for i in 1:n
    for j in i:n
        # Extract columns as vectors
        x = td_by_simrep_by_region[!, i]
        y = td_by_simrep_by_region[!, j]
        # Find indices where both are not missing
        #idx = .!ismissing.(x) .& .!ismissing.(y)
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
        corrmat_idx[i, j] = sum(idx)
        corrmat_idx[j, i] = sum(idx) # Symmetric
    end
end
#println(corrmat_idx)
#mean(corrmat_idx) # 31.0 for psampled = 0.1

corr_df = DataFrame(corrmat, Symbol.(names(td_by_simrep_by_region)))
#rename!(corr_df, Symbol.(names(td_by_simrep_by_region)))
maximum(skipmissing(corrmat))
maximum(x for x in skipmissing(corrmat) if x != 1.0)
minimum(skipmissing(corrmat))
# CHECK - for psampled = 1.0, minimum correlation is -0.52 between TLK3 (Cornwall & Isles of Scilly)
# and TLK5 (West of England)
# Only includes TDs for simulations in which both regions had at least 1 ICU case 
idx = isfinite.(td_by_simrep_by_region.TLK3) .& isfinite.(td_by_simrep_by_region.TLK5)
scatter( td_by_simrep_by_region.TLK3[idx]
    , td_by_simrep_by_region.TLK5[idx]
    , xlabel = "TLK3 TD (days)", ylabel = "TLK5 TD (days)", label = "")
plot( td_by_simrep_by_region.TLK3[idx], xlabel = "simulation replicate", ylabel = "TD (days)", label = "TLK3: Cornwall & Isles of Scilly" )
plot!( td_by_simrep_by_region.TLK5[idx] , label = "TLK5: West of England" )
# CHECK - for psampled = 0.1, maximum correlation is 1.00 including between TLD1 (Cumbria)
# and TLK4 (Devon)
# Only includes TDs for simulations in which both regions had at least 1 ICU case 
idx = isfinite.(td_by_simrep_by_region.TLD1) .& isfinite.(td_by_simrep_by_region.TLK4)
scatter( td_by_simrep_by_region.TLD1[idx]
    , td_by_simrep_by_region.TLK4[idx]
    , xlabel = "TLD1:Cumbria TD (days)", ylabel = "TLK4:Devon TD (days)", label = "")
plot( td_by_simrep_by_region.TLD1[idx], xlabel = "simulation replicate", ylabel = "TD (days)", label = "TLD1: Cumbria" )
plot!( td_by_simrep_by_region.TLK4[idx] , label = "TLK4: Devon" )


# Plot correlations
# Create region codes/names for labelling heatmap
x_names = names(td_by_simrep_by_region)
ITL2_dict = Dict(zip(REGKEY.code, REGKEY.name))
ITL2_names = [ITL2_dict[k] for k in x_names]
# Merge pairs into strings (e.g., "key:value" = "code:names")
y_names = [string(k, ":", v) for (k, v) in zip(x_names, ITL2_names)]
#println(y_names)

# Create correlation values for annotating the heatmap (if required)
anns = []
for i in 1:n
    for j in 1:n
        val = round(corrmat[i, j], digits=2)
        # Note: y-axis is reversed in Plots.jl heatmap to match matrix layout
        push!(anns, (j, n - i + 1, string(val)))
    end
end
# Plot heatmap
heatmap(  x_names # x-tick labels (columns)
        , xmirror = true # move x-asix labels to the top
        , xrotation = 90 # Rotate x-axis labels
        #, xticks = (1:n, names(td_by_simrep_by_region))
        , xticks = (collect(1:n) .- 0.5, x_names) #names(td_by_simrep_by_region)) # centers of each cell
        , y_names #names(td_by_simrep_by_region) # y-tick labels (rows)
        #, yticks = (1:n, names(td_by_simrep_by_region))
        , yticks = (collect(1:n) .- 0.5, y_names) #names(td_by_simrep_by_region))# centers of each cell
        , corrmat;              # the correlation matrix
          color = :RdBu        # diverging color map, blue to red
        , clim = (-1, 1)       # color limits for correlation coefficients
        , xlabel = "ITL2 region"
        , ylabel = "ITL2 region"
        , title = "Correlation Matrix for regional TD across 1000 simulation replicates"
        , yflip = true         # flips y-axis to match matrix layout
        , size = (1100, 800)
        #, annotations = anns
)

# Save correlation matrix to a file
corr_df_w_names = insertcols!(corr_df, 1, :ITL2_region => y_names )
println(corr_df_w_names)
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/regional_TD_correlation_matrix_p_0.1.csv", corr_df_w_names)
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/regional_TD_correlation_matrix_p_1.csv", corr_df_w_names)
# for psampled = 0.1
CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/psampled_10/regional_TD_correlation_matrix_p_10.csv", corr_df_w_names)


#println(countmap(corr_df_w_names[:,2:size(corr_df_w_names,2)])) 
#sort( countmap( corrmat ) )

# List region pairs with a correlation above a specified threshold
#corr_matrix = cor(Matrix(df))
varnames = names(td_by_simrep_by_region)
n = length(varnames)
if psampled == 0.1
    threshold = 0.75
elseif psampled == 1.0
    threshold = 0.25
end

results = DataFrame(Row = String[], Col = String[], Correlation = Float64[])

for i in 1:n-1
    for j in i+1:n
        val = corrmat[i, j]
        if ismissing(val)
            continue
        else
            if abs(val) > threshold
                push!(results, (varnames[i], varnames[j], val))
            end
        end
    end
end

sort!(results, :Correlation, rev=true)#, by=abs)
println(results)

# Build a Dict for fast lookup
region_lookup = Dict(REGKEY.code .=> REGKEY.name)

# Prepare results DataFrame with extra columns
results = DataFrame(Region_1 = String[], Region_2 = String[]
                    , Region_1_name = String[], Region_2_name = String[]
                    , Correlation = Float64[]
                    , sim_replicates_n = Int[]
                    #, Stat_Power = Float64[], required_n = Int[] 
                    )

for i in 1:n-1
    for j in i+1:n
        val = round( corrmat[i, j], digits = 2)
        if ismissing(val)
            continue
        else
            if abs(val) > threshold
                row_key = varnames[i]
                col_key = varnames[j]
                row_val = get(region_lookup, row_key, missing)
                col_val = get(region_lookup, col_key, missing)
                #stat_power = round( corr_power_region_pair[i,j], digits = 2) #TODO link this to the actual correlation value between region pairs rather than set value for all pairs
                sim_replicates_n = Int( ceil( corrmat_idx[i,j] ) )
                #req_n = Int( ceil( required_n_region_pair[i,j] ) ) #TODO link this to the actual correlation value between region pairs rather than set value for all pairs
                push!(results, (row_key, col_key, row_val, col_val, val
                                , sim_replicates_n
                                #, stat_power, , req_n 
                                )
                     )
            end
        end
    end
end

# Sort results by absolute correlation, high to low
sort!(results, :Correlation, rev=true)#, by=abs)
println(results)

#TEST
# Frequency table of homeregions
#println(countmap(G.homeregion))


### For each region, fit gamma distribution to infection times (tinf) corresponding the minimum report time (treport / TD)

# Function to compute the negative log-likelihood of a gamma distribution.
# This can then be optimised to fit truncated data (the tinf times which are truncated at the simulation maxtime).

#params = init_gamma

# Negative log-likelihood for truncated Gamma
#function nll_trunc_gamma(params)
function nll_trunc_gamma(; params, times, trunc)
    shape, scale = params
    if shape <= 0 || scale <= 0
        return 1e10 #Inf
    end
    
    # Define distribution
    d0 = Gamma(shape, scale)
    #norm = cdf(d0, upper) # Only upper truncation
    #ll = sum(logpdf(d0, t) for t in times) - length(times) * log(norm)
    
    # Define truncated distribution
    lower, upper = trunc
    d = truncated( d0, lower, upper)
    
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    # Compute log-likelihood
    ll = sum( sort(logpdf.(d, times_no_zeros)) ) #ll = sum( sort(logpdf.(d, times)) ) 
          
    return isfinite(-ll) ? -ll : 1e10 #return -ll
end


# Negative log-likelihood for truncated Weibull
#function nll_trunc_weibull(params)
function nll_trunc_weibull(; params, times, trunc)
    α, θ = params
    if α <= 0 || θ <= 0
        return 1e10 #Inf
    end
    
    # Define distribution
    d0 = Weibull(α, θ)
    
    # Define truncated distribution
    lower, upper = trunc
    d = truncated(d0, lower, upper)
    
    #if any(x -> x < lower || x > upper, times)
    #    return 1e10 #Inf
    #end
    
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    
    # Compute log-likelihood
    ll = sum(logpdf.(d, times_no_zeros)) #ll = sum(logpdf.(d, times))
    
    return isfinite(-ll) ? -ll : 1e10 #return -ll
end

#histogram(times)

# TEST function to plot data and fitted distribution
function dist_fit_plot(; d0 = Gamma( res_gamma.minimizer[1], res_gamma.minimizer[2])
                        , lower = 0.0, upper = 60.0
                        , data_to_fit = times
                        , data_type = "tinf" #"(treport-tinf)
                        )
    # Convert distribution to truncated version
    d = truncated( d0, lower, upper)
    # sample data points from fitted truncated distribution
    #d0_sample=rand(d0,10000)
    d_sample=rand(d,1000000)
        
    # Plot histogram of sampled data
    histogram( d_sample, alpha=0.5, normalize=:pdf, label="samples drawn from fitted dist")
    #histogram!(d0_sample,alpha=0.5, normalize=:pdf, color=:yellow)
    
    # Add simulated outbreak data to histogram
    histogram!(data_to_fit, normalize=:pdf, color=:red, alpha=0.5, label=data_type*" data")#" tinf or (treport-tinf) data")
    
    # Add simulated outbreak data with zero values removed to histogram
    data_to_fit_no_zeros = sort( filter(x -> x != 0, data_to_fit) )
    histogram!(data_to_fit_no_zeros, normalize=:pdf, color=:green, alpha=0.5, label=data_type*" data")#" tinf or (treport-tinf) data")
    
    # Overlay fitted (truncated) distribution PDF
    xs = range(0, upper; length=300)
    norm = cdf(d, upper) # normalization for truncation
    t_pdf_vals = pdf.(d, xs) ./ norm # truncated PDF
    plot!(xs, t_pdf_vals; label=data_type*" truncated fit", color=:blue, lw=2)
    
end

# Function to compare distribution fits
function fit_multi_dist(; times, init_gamma, init_nbinom, init_lognorm, init_weibull, lower, upper)
    
    # Fit all distributions and record results in a dictionary
    fit_results = Dict()

    # Gamma
    obj_f_gamma = init_gamma -> nll_trunc_gamma(;params = init_gamma, times = times, trunc = [0, 60])
    #init_gamma=[16.466863607650694,2.8170337707969924]
    #res_gamma = optimize(nll_trunc_gamma, init_gamma, NelderMead() #GoldenSection()) #Brent()) #AcceleratedGradientDescent()) #MomentumGradientDescent()) #GradientDescent()) #ConjugateGradient()) #LBFGS()) #SimulatedAnnealing()) #NelderMead())
    res_gamma = optimize(obj_f_gamma, init_gamma, NelderMead() #GoldenSection()) #Brent()) #AcceleratedGradientDescent()) #MomentumGradientDescent()) #GradientDescent()) #ConjugateGradient()) #LBFGS()) #SimulatedAnnealing()) #NelderMead())
                        , Optim.Options( iterations = 2000    # maximum number of iterations
                                         , g_tol = 1e-8         # gradient tolerance
                                         , store_trace = true   # store the optimization trace
                                         #, show_trace = true    # print progress to stdout
                                         , show_warnings = true  # show warnings
                                        )
                        ) 
    
    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converged( res_gamma )
        #nll_g = nll_trunc_gamma( Optim.minimizer( res_gamma ) )
        nll_g = nll_trunc_gamma( params = Optim.minimizer( res_gamma ), times = times, trunc = [lower, upper] )
        k = length( Optim.minimizer( res_gamma ) ) # Number of parameters estimated
        aic_g = 2*k + 2*nll_g 
        fit_results["Gamma"] = ( nll_g, aic_g, Optim.minimizer( res_gamma ) )
    #else
    #    fit_results["Gamma"] = ("NA", "NA", "NA")
    #end

    # TEST PLOT TO VIEW DATA AND FIT
    # dist_fit_plot(; d0 = Gamma( res_gamma.minimizer[1], res_gamma.minimizer[2]), lower = 0.0, upper = 60.0, data_to_fit = times, data_type = "tinf")
    
    # Weibull
    obj_f_weibull = init_weibull -> nll_trunc_weibull(;params = init_weibull, times = times, trunc = [lower, upper])
    #res_weibull = optimize(nll_trunc_weibull, init_weibull, NelderMead()
    res_weibull = optimize(obj_f_weibull, init_weibull, NelderMead() 
                            , Optim.Options(
                                              iterations = 2000    # maximum number of iterations
                                            , g_tol = 1e-8         # gradient tolerance
                                            , store_trace = true   # store the optimization trace
                                            #, show_trace = true    # print progress to stdout
                                            , show_warnings = true  # show warnings
                                            )
    )
    
    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converge( res_weibull ) #res_weibull.exitflag == :Success
        #nll_w = nll_trunc_weibull(Optim.minimizer(res_weibull))
        nll_w = nll_trunc_weibull( params = Optim.minimizer( res_weibull ), times = times, trunc = [lower, upper] )
        k = length(Optim.minimizer(res_weibull)) # Number of parameters estimated
        aic_w = 2*k + 2*nll_w
        fit_results["Weibull"] = (nll_w, aic_w, Optim.minimizer(res_weibull))
    #else
    #   fit_results["Weibull"] = ("NA", "NA", "NA")
    #end

    # TEST PLOT TO VIEW DATA AND FIT
    # dist_fit_plot(; d0 = Weibull( res_weibull.minimizer[1], res_weibull.minimizer[2]), lower = 0.0, upper = 60.0, data_to_fit = times, data_type = "tinf")

    # Report results
    if isempty(fit_results)
        println("All fits failed.")
    else
        println("Goodness-of-fit (AIC) for each distribution:")
        for (dist, (_nll, aic, params)) in fit_results
            println("  $dist: negative log-likelihood = $_nll, AIC = $aic, parameters = $params")
        end
    # Select best distribution
    best_dist = findmin([(aic, dist) for (dist, (_nll, aic, _params)) in fit_results])[1][2]
    println("Best-fitting distribution: $best_dist ") #$fit_results[$("$best_dist")]")
    end

    return( [fit_results] )

end

# Simulation maxtime
maxtime = 60.0

# Initialise df to record regions and the median tinf and treport (TD) times
median_times_by_region = DataFrame( 
      ITL2_code = String[]
    , ITL2_name = String[]
    , tinf_median = Float64[]
    , tinf_median_fit = Float64[]
    , treport_median = Float64[]
    , treport_median_fit = Float64[]
)
# Loop through all ITL2 regions

for r in REGKEY.code # TEST: r = "TLI7"; r = icu_regions[2]

    # Times to fit truncated gamma distribution to
    times = tinf_by_region_dict[ r ]
    #Test
    #histogram(times)
    # CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/tinf_times_TLI7.csv", DataFrame(value = times))

    # Define truncation range
    lower = 0 #0.00001
    upper = maxtime

    # Initial guesses
    #init_gamma = [10.0, 5.0] # shape*scale = mean. Median is ~50 for the tinf data # [1.0, 1.0]
    init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
    #init_gamma_mixed = [init_gamma; 0.05] # For use in a gamma mixed model with probability of zero
    init_weibull = init_gamma #[20.0, 80.0]
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    
    fit_tinf_results = fit_multi_dist(; times
                                    , init_gamma, init_nbinom, init_lognorm, init_weibull
                                    , lower, upper)
    
    # Convert to a df and print
    fit_tinf_results_df_rows = [
                                (
                                  dist
                                , v[1] 
                                , v[2]
                                , v[3][1]
                                , v[3][2]
                                , missing  # spare column, can be changed as needed
                                ) for (dist, v) in fit_tinf_results[1]
                                ]
    fit_tinf_results_df = DataFrame(fit_tinf_results_df_rows
                                    , [:Distribution, :Negative_log_lik, :AIC, :param1, :param2, :param3])
    println(fit_tinf_results_df)

    ## Compute fit data for plotting
    # Gamma
    gamma_shape_fit, gamma_scale_fit = fit_tinf_results[1]["Gamma"][3]
    tinf_fitted_gamma = Gamma(gamma_shape_fit, gamma_scale_fit)
    tinf_fitted_gamma_trunc = truncated(tinf_fitted_gamma, 0.0, upper)
    # Overlay fitted (truncated) gamma PDF
    xs = range(0, upper; length=300)
    norm_gamma = cdf(tinf_fitted_gamma, upper) # normalization for truncation
    tinf_pdf_vals_gamma = pdf.(tinf_fitted_gamma, xs) ./ norm_gamma # truncated PDF
    tinf_pdf_vals_gamma_trunc = pdf.(tinf_fitted_gamma_trunc, xs) #./ norm_gamma # truncated PDF
    #CHECK
    #plot(tinf_pdf_vals_gamma)
    #plot!(tinf_pdf_vals_gamma_trunc, col="red")
    
    # Weibull
    weibull_shape_fit, weibull_scale_fit = fit_tinf_results[1]["Weibull"][3]
    tinf_fitted_weibull = Gamma(weibull_shape_fit, weibull_scale_fit)
    tinf_fitted_weibull_trunc = truncated(tinf_fitted_weibull, 0.0, upper)
    # Overlay fitted (truncated) gamma PDF
    norm_weibull = cdf(tinf_fitted_weibull, upper) # normalization for truncation
    tinf_pdf_vals_weibull = pdf.(tinf_fitted_weibull, xs) ./ norm_weibull # truncated PDF
    tinf_pdf_vals_weibull_trunc = pdf.(tinf_fitted_weibull_trunc, xs) #./ norm_weibull # truncated PDF
    #CHECK
    #plot(tinf_pdf_vals_weibull)
    #plot!(tinf_pdf_vals_weibull_trunc, col="red")
    
    # Plot histogram (normalized to probability density)
    p1 = histogram(times; bins=30
                    , normalize=:pdf
                    , label="tinf"
                    , color=:lightblue, alpha=0.5, legend=:topright
                    , title="$(r) $(REGKEY[REGKEY.code .== r, :name][1])" 
                    , xlabel="Days", ylabel = "Density"
                    , xlims = [0,150]
                    , ylims = [0,0.15]
                   )

    # Plot tinf fit
    p1 = plot!(xs, tinf_pdf_vals_gamma; label="tinf truncated Gamma fit", color=:blue, lw=2)
    p1 = plot!(xs, tinf_pdf_vals_weibull; label="tinf truncated Weibull fit", color=:purple, lw=2)

    # Add median for tinf
    tinf_median = quantile(times, 0.5)
    times_no_zeros = sort( filter(x -> x != 0, times) )
    tinf_median_no_zeros = quantile(times_no_zeros, 0.5)
    #tinf_median = median(times)
    p1 = vline!([tinf_median], color=:blue
            , label="tinf median")
    p1 = vline!([tinf_median_no_zeros], color=:purple
            , label="tinf (no zeros) median")

    # Add median for fitted tinf
    tinf_median_fit_gamma = quantile(tinf_fitted_gamma_trunc, 0.5)
    p1 = vline!([tinf_median_fit_gamma], color=:blue, linestyle=:dash
            , label="tinf median truncated Gamma fit")
    tinf_median_fit_weibull = quantile(tinf_fitted_weibull_trunc, 0.5)
    p1 = vline!([tinf_median_fit_weibull], color=:purple, linestyle=:dash
            , label="tinf median truncated Weibull fit")

    # Add histogram for treport
    p1 = histogram!(treport_by_region_dict[ r ]; bins=30
                , normalize=true
                , label="treport", color=:pink, alpha=0.5)

    # Add median for treport
    treport_median = quantile(treport_by_region_dict[ r ], 0.5)
    p1 = vline!([treport_median], color=:red, label="treport median")
    
    # Fit gamma distribution to (treport - tinf) from earliest treport
    times = tinf_treport_by_region_dict[ r ]
    
    # Initial guesses
    #init_params = [1.0, 1.0]
    init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
    init_weibull = init_gamma
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    
    fit_tinf_treport_results = fit_multi_dist(; times
                                               , init_gamma, init_nbinom, init_lognorm, init_weibull
                                               , lower, upper)
    println(fit_tinf_treport_results)

    # Gamma fit
    tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit = fit_tinf_treport_results[1]["Gamma"][3]
    tinf_treport_fitted_gamma = Gamma(tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit)
    tinf_treport_median_gamma_fit = median( tinf_treport_fitted_gamma ) # Note that this distribution is not truncated
    
    # Add median for fitted truncated gamma tinf + fitted gamma (treport-tinf) = fitted treport
    p1 = vline!([tinf_treport_median_gamma_fit + tinf_median_fit_gamma], color=:red, linestyle=:dash
            , label="treport median Gamma fit")

    # Add histogram for (treport - tinf)
    p2 = histogram(tinf_treport_by_region_dict[ r ]; bins=30, normalize=:pdf
                   , label="treport-tinf", color=:green, alpha=0.5, legend=:topright
                   , xlabel="Days", ylabel="Density"
                   , xlims = [0,150]
                   , ylims = [0,0.15]
                    )

    
    # Overlay fitted gamma PDF
    tinf_treport_pdf_vals_gamma = pdf.(tinf_treport_fitted_gamma, xs)
    # Plot (treport - tinf) fit
    p2 = plot!(xs, tinf_treport_pdf_vals_gamma; label="treport-tinf gamma fit"
            , color=:green, lw=2)

    # Add median for (treport - tinf)
    tinf_treport_median = quantile(tinf_treport_by_region_dict[ r ], 0.5) # Note does not include any zeros
    p2 = vline!([tinf_treport_median], color=:green
            , label="treport-tinf median", lw = 2)

    # Add median for fitted (treport - tinf)
    tinf_treport_median_fit_gamma = quantile(tinf_treport_fitted_gamma, 0.5)
    p2 = vline!([tinf_treport_median_fit_gamma], color=:green, linestyle=:dash
            , label="treport-tinf median gamma fit")

    plot( size=(800, 1200), p1, p2, layout = (2,1))
    
    #savefig("scripts/median_TD_by_region_plots/$r.png")
    #savefig("scripts/median_TD_by_region_plots/sim_regionentry/$r.png")
    # For psampled = 1
    #savefig("scripts/median_TD_by_region_plots/sim_regionentry/psampled_100/$r.png")
    # For psampled = 1
    savefig("scripts/median_TD_by_region_plots/sim_regionentry/psampled_10/$r.png")

    # Add row to df to record the median value by region
    new_row = ( ITL2_code = r
                , ITL2_name = REGKEY[REGKEY.code .== r, :name][1]
                , tinf_median = tinf_median
                , tinf_median_fit = tinf_median_fit_gamma
                , treport_median = treport_median
                , treport_median_fit = (tinf_treport_median_fit_gamma + tinf_median_fit_gamma)
                )
    push!( median_times_by_region, new_row )
     
end

println( median_times_by_region )
median_times_by_region_sorted = sort( median_times_by_region , :treport_median_fit )
println( median_times_by_region_sorted )
#CSV.write("scripts/median_TD_by_region_plots/TD_by_region_sorted.csv", median_times_by_region_sorted)
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/TD_by_region_sorted.csv", median_times_by_region_sorted)
# For psampled = 1
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/psampled_100/TD_by_region_sorted_psampled_100.csv", median_times_by_region_sorted)
# For psampled = 0.1
CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/psampled_10/TD_by_region_sorted_psampled_10.csv", median_times_by_region_sorted)

# Add the percent of simulation replicates that include ICU cases for each region
median_times_by_region_sorted[!,:perc_simreps_have_icu_cases] = zeros(nrow(median_times_by_region_sorted))
for i in 1:nrow( median_times_by_region_sorted )
    median_times_by_region_sorted.perc_simreps_have_icu_cases[i] = perc_simreps_icu_by_region[ 1 , median_times_by_region_sorted[i,"ITL2_code"] ]
end    
println(median_times_by_region_sorted)
minimum(median_times_by_region_sorted.perc_simreps_have_icu_cases )
maximum(median_times_by_region_sorted.perc_simreps_have_icu_cases )

# For each region find the top 3 regions with the highest correlations and add name and correlation value
# Add columns
for i in 1:3
    median_times_by_region_sorted[!, Symbol("high_cor_$(i)_region")] = fill("", nrow(median_times_by_region_sorted))
    median_times_by_region_sorted[!, Symbol("high_cor_$(i)_correlation")] = zeros(nrow(median_times_by_region_sorted))
end
# Fill columns
for i in 1:length(median_times_by_region_sorted.ITL2_code)
    ## Find regions with the top 3 highest correlations with region of interest
    # Region correlations with region i
    region_correls = DataFrame( ITL2_region = corr_df_w_names[:,1]
                                , correlation = corr_df_w_names[:,median_times_by_region_sorted.ITL2_code[i]]
                              )
    # Remove correlation with self
    region_i = median_times_by_region_sorted.ITL2_code[i]
    self_index = findfirst(j -> first(region_correls.ITL2_region[j], 4) == first(region_i, 4)
                            , 1:nrow(region_correls))
    # Remove self-index row
    if self_index !== nothing
        region_correls_ex_self = filter(row -> row !== self_index, eachindex(eachrow(region_correls))) .|> x -> region_correls[x, :]
        region_correls_ex_self = DataFrame(region_correls_ex_self)  # Reconstruct df2 from filtered rows
        #println(region_correls_ex_self)
    end
    #println(region_correls_ex_self.ITL2_region )
    # Filter out regions with missing correlation (i.e. no overlapping simulation replicates with region i)
    region_correls_filtered = region_correls_ex_self[ findall(!ismissing, region_correls_ex_self[:,2] ), : ]
    # Select top 3 correlations with region i
    top_3 = sort( region_correls_filtered , :correlation , rev=true)[1:3,1:2]
    # Fill data in median_times_by_region_sorted df
    for k in 1:3
        median_times_by_region_sorted[!, Symbol("high_cor_$(k)_region")][i]      = top_3[k,1]
        median_times_by_region_sorted[!, Symbol("high_cor_$(k)_correlation")][i] = top_3[k,2]
    end
       
end

println(median_times_by_region_sorted)

#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/TD_by_region_sorted.csv", median_times_by_region_sorted)
# For psampled = 1 #=100%
#CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/psampled_100/TD_by_region_sorted_psampled_100.csv", median_times_by_region_sorted)
# For psampled = 0.1 #=10%
CSV.write("scripts/median_TD_by_region_plots/sim_regionentry/psampled_10/TD_by_region_sorted_psampled_10.csv", median_times_by_region_sorted)



# TODO possibly resimulate
### Apply fitted statistical distribution to the minimum report time (treport / TD)

# Resimulate tsample and treport based on truncated fit to tinf

# Draw 1000 tinf from fitted distribution
tinf_fit_samples = rand(tinf_trunc_gamma, 1000)
histogram!(tinf_fit_samples, alpha = 0.1)
histogram!(tinf_fit_samples; bins=30, normalize=true, label="tinf fit samples", color=:lightgreen, alpha=0.6, legend=:topleft)

# Re-simulate tsample and treport


