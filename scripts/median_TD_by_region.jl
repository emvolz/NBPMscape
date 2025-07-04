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
#using StatsBase
using GLM, Statistics
using Distributions
using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using DataFrames
import UUIDs 
import StatsBase 
using Interpolations
import SpecialFunctions as SF 
using Plots 
using LinearAlgebra
using Optim

using RData 
using CSV 

# Load simulation
sims = load("covidlike-1.0-sims.jld2" , "sims")

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
psampled = 1
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

# Loop through each simulation
for s in 1:length(sims)
    #Test
    #s=1
    
    fo = sims[s]
    
    # Filter for ICU cases 
    G = fo.G[ isfinite.(fo.G.ticu), : ]
    
    # List regions with an ICU admission (i.e. possibility of getting a sample and detection)
    icu_regions = unique(G.homeregion)
    
    # Loop through each region and record time of infection stats
    for r in icu_regions #REGKEY.code
        #Test
        #r=icu_regions[1]
        
        g_region = G[ (G.homeregion).==r, :]
        
        # Sample size 
	    n = rand( Binomial( size(g_region,1), psampled) ) #p.psampled ) )
	    # Subforest within region
	    g_region_sub = g_region[sample( 1:size(g_region,1), n, replace=false ), :]

        # Check if there are any rows before continuing
        if size(g_region_sub,1) == 0 
            continue
        elseif size(g_region_sub,1) > 0 
            
            #= Sample time has uniform distribution between time of admission to ICU and time of recovery
               TODO THIS MAY NEED TO BE UPDATED FOR LATER VERSIONS WITH MORE COMPLEX CARE PATHWAYS
               NOTE TIME BETWEEN ticu and trecovered can be large and so with uniform distribution of 
               sampling the tsample can be a long time after tinf. Example I saw was tinf = 48.9, ticu = 52.6, trecovered = 91.1, and treport = 90.7
               If we're modelling 100% of ICU cases being sampled then more likely to sampled closer to ticu
            =#
            
            # Generate sample times
            #tsample = (size(g_region,1)>0) ? map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region) ) : []
            tsample = map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region_sub) )
            g_region_sub.tsample = tsample 
        
            # Simulate reports times
            #treport = (size(g_region,1)>0) ? (tsample .+ turnaroundtime) : []
            treport = (tsample .+ turnaroundtime)
            g_region_sub.treport = treport
        
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
            
        end
                
    end

end

#TEST
# Frequency table of homeregions
#println(countmap(G.homeregion))


### For each region, fit gamma distribution to infection times (tinf) corresponding the minimum report time (treport / TD)

# Function to compute the negative log-likelihood of a gamma distribution.
# This can then be optimised to fit truncated data (the tinf times which are truncated at the simulation maxtime).
function negloglik_gamma(params)
    shape, scale = params
    if shape <= 0 || scale <= 0
        return Inf
    end
    d0 = Gamma(shape, scale)

    #norm = cdf(d0, upper) # Only upper truncation
    #ll = sum(logpdf(d0, t) for t in times) - length(times) * log(norm)

    d = truncated( d0, 0, upper)
    ll = sum(logpdf.(d, times))

    return -ll
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

for r in REGKEY.code #icu_regions
#TEST
#r = "TLE3"
#r = icu_regions[2]

    # Times to fit truncated gamma distribution to
    times = tinf_by_region_dict[ r ]
    #Test
    histogram(times)

    # Define truncation time
    upper = maxtime

    # Fit truncated gamma distribution to tinf from earliest treport
    init_params = [1.0, 1.0]
    result = optimize(negloglik_gamma, init_params, BFGS())
    shape_hat, scale_hat = Optim.minimizer(result)

    tinf_fitted_gamma = Gamma(shape_hat, scale_hat)
    tinf_trunc_gamma = truncated(tinf_fitted_gamma, 0.0, upper)

    # Plot histogram (normalized to probability density)
    p1 = histogram(times; bins=30, normalize=:pdf, label="tinf"
             , color=:lightblue, alpha=0.5, legend=:topright
             , title="$(r) $(REGKEY[REGKEY.code .== r, :name][1])" 
             , xlabel="Days", ylabel = "Density"
             )

    # Overlay fitted (truncated) gamma PDF
    xs = range(0, upper; length=300)
    norm = cdf(tinf_fitted_gamma, upper) # normalization for truncation
    tinf_pdf_vals = pdf.(tinf_fitted_gamma, xs) ./ norm # truncated PDF

    # Plot tinf fit
    p1 = plot!(xs, tinf_pdf_vals; label="tinf truncated gamma fit", color=:blue, lw=2)

    # Add median for tinf
    tinf_median = quantile(times, 0.5)
    p1 = vline!([tinf_median], color=:blue
            , label="tinf median")

    # Add median for fitted tinf
    tinf_median_fit = quantile(tinf_trunc_gamma, 0.5)
    p1 = vline!([tinf_median_fit], color=:blue, linestyle=:dash
            , label="tinf median truncated gamma fit")

    # Add histogram for treport
    p1 = histogram!(treport_by_region_dict[ r ]; bins=30, normalize=true
                , label="treport", color=:pink, alpha=0.5)

    # Add median for treport
    treport_median = quantile(treport_by_region_dict[ r ], 0.5)
    p1 = vline!([treport_median], color=:red, label="treport median")
    
    # Fit gamma distribution to (treport - tinf) from earliest treport
    times = tinf_treport_by_region_dict[ r ]
    init_params = [1.0, 1.0]
    result = optimize(negloglik_gamma, init_params, BFGS())
    shape_hat, scale_hat = Optim.minimizer(result)

    tinf_treport_fitted_gamma = Gamma(shape_hat, scale_hat)
    
    # Add median for fitted truncated gamma tinf + fitted gamma (treport-tinf) = fitted treport
    p1 = vline!([tinf_treport_median_fit + tinf_median_fit], color=:red, linestyle=:dash
            , label="treport median fit")

    # Add histogram for (treport - tinf)
    p2 = histogram(tinf_treport_by_region_dict[ r ]; bins=30, normalize=:pdf
                   , label="treport-tinf", color=:green, alpha=0.5, legend=:topright
                   , xlabel="Days", ylabel="Density")

    # Overlay fitted gamma PDF
    xs = range(0, upper; length=300)
    #norm = cdf(tinf_fitted_gamma, upper) # normalization for truncation
    tinf_treport_pdf_vals = pdf.(tinf_treport_fitted_gamma, xs) #./ norm # truncated PDF

    # Plot (treport - tinf) fit
    p2 = plot!(xs, tinf_treport_pdf_vals; label="treport-tinf gamma fit"
            , color=:green, lw=2)

    # Add median for (treport - tinf)
    tinf_treport_median = quantile(tinf_treport_by_region_dict[ r ], 0.5)
    p2 = vline!([tinf_treport_median], color=:green
            , label="treport-tinf median", lw = 2)

    # Add median for fitted (treport - tinf)
    tinf_treport_median_fit = quantile(tinf_treport_fitted_gamma, 0.5)
    p2 = vline!([tinf_treport_median_fit], color=:green, linestyle=:dash
            , label="treport-tinf median gamma fit")

    plot( size=(800, 1200), p1, p2, layout = (2,1))
    
    savefig("scripts/median_TD_by_region_plots/$r.png")

    # Add row to df to record the median value by region
    new_row = ( ITL2_code = r
                , ITL2_name = REGKEY[REGKEY.code .== r, :name][1]
                , tinf_median = tinf_median
                , tinf_median_fit = tinf_median_fit
                , treport_median = treport_median
                , treport_median_fit = (tinf_treport_median_fit + tinf_median_fit)
                )
    push!( median_times_by_region, new_row )
     
end

println( median_times_by_region )
println( sort( median_times_by_region , :treport_median_fit )  )

# TODO possibly resimulate
### Apply this fitted gamma distribution to the minimum report time (treport / TD)

# Resimulate tsample and treport based on trunacated gamma fit to tinf

# Draw 1000 tinf from fitted distribution
tinf_fit_samples = rand(tinf_trunc_gamma, 1000)
histogram!(tinf_fit_samples, alpha = 0.1)
histogram!(tinf_fit_samples; bins=30, normalize=true, label="tinf fit samples", color=:lightgreen, alpha=0.6, legend=:topleft)

# Re-simulate tsample and treport


