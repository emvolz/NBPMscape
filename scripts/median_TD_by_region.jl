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
using Interpolations
import SpecialFunctions as SF 
using Plots 
using LinearAlgebra
using Optim
using RData 
using CSV 

# Load simulation
sims = load("covidlike-1.0-sims.jld2" , "sims") # preliminary simulation
#sims = load("covidlike-1.0-sims-regionentry.jld2" , "sims") # preliminary simulation re-run with regionentry adjusted

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

function nll_trunc_gamma(params)
    shape, scale = params
    if shape <= 0 || scale <= 0
        return 1e10 #Inf
    end
    d0 = Gamma(shape, scale)
    #norm = cdf(d0, upper) # Only upper truncation
    #ll = sum(logpdf(d0, t) for t in times) - length(times) * log(norm)
    d = truncated( d0, lower, upper)
    #if any(x -> x < lower || x > upper, sort(times))
    #    return 1e100 #Inf
    #end

    # Filter for finite values
    #logpdf_vec = sort(logpdf.(d, times))
    #logpdf_finite_vec = filter(isfinite, logpdf_vec)
    #ll = sum( logpdf_finite_vec )
    
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    # Compute log-likelihood
    ll = sum( sort(logpdf.(d, times_no_zeros)) ) #ll = sum( sort(logpdf.(d, times)) ) 
          
    return isfinite(-ll) ? -ll : 1e10 #return -ll
end

#TEST
#test_d=rand(d,10000)
#histogram(test_d,alpha=0.5, normalize=:pdf)
#test_d0=rand(d0,10000)
    #histogram!(test_d0,alpha=0.5, normalize=:pdf, color=:green)
    #histogram!(times, normalize=:pdf, color=:red, alpha=0.5)
    
    # Overlay fitted (truncated) gamma PDF
    #xs = range(0, upper; length=300)
    #norm = cdf(d, upper) # normalization for truncation
    #tinf_pdf_vals = pdf.(d, xs) ./ norm # truncated PDF

    # Plot tinf fit
    #plot!(xs, tinf_pdf_vals; label="tinf truncated gamma fit", color=:blue, lw=2)
    #plot!(xs, pdf.(Gamma(18,3),xs); label="tinf truncated gamma fit", color=:blue, lw=2)


function nll_trunc_nbinom(params)
    r, p = params
    if r <= 0 || p <= 0
        return 1e10 #Inf
    end
    d0 = NegativeBinomial(r, p)

    d = truncated( d0, lower, upper)
    
    #if any(x -> x < lower || x > upper, sort(times))
    #    return 1e100 #Inf
    #end
    
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    # Compute log-likelihood
    #ll = sum( sort(logpdf.(d, times)) )
    ll = sum( sort(logpdf.(d, times_no_zeros)) ) 
    
    return isfinite(-ll) ? -ll : 1e10
end

# Negative log-likelihood for truncated log-normal
function nll_trunc_lognorm(params)
    μ, σ = params
    lower = 0
    if σ <= 0
        return 1e10 #Inf
    end
    d0 = LogNormal(μ, σ)
    d = truncated(d0, lower, upper)
    if any(x -> x < lower || x > upper, times)
        return 1e10 #Inf
    end
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    # Compute log-likelihood
    ll = sum(logpdf.(d, times_no_zeros)) #ll = sum(logpdf.(d, times))
    return isfinite(-ll) ? -ll : 1e10
end

# Negative log-likelihood for truncated Weibull
function nll_trunc_weibull(params)
    α, θ = params
    lower = 0
    if α <= 0 || θ <= 0
        return 1e10 #Inf
    end
    d0 = Weibull(α, θ)
    d = truncated(d0, lower, upper)
    if any(x -> x < lower || x > upper, times)
        return 1e10 #Inf
    end
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    # Compute log-likelihood
    ll = sum(logpdf.(d, times_no_zeros)) #ll = sum(logpdf.(d, times))
    return isfinite(-ll) ? -ll : 1e10 #return -ll
end

histogram(times)

# Function to compare distribution fits
function fit_multi_dist(; times, init_gamma, init_nbinom, init_lognorm, init_weibull, lower, upper)
    
    # Fit all distributions and record results in a dictionary
    fit_results = Dict()

    # Gamma
    res_gamma = optimize(nll_trunc_gamma, init_gamma, NelderMead() #GoldenSection()) #Brent()) #AcceleratedGradientDescent()) #MomentumGradientDescent()) #GradientDescent()) #ConjugateGradient()) #LBFGS()) #SimulatedAnnealing()) #NelderMead())
                        , Optim.Options( iterations = 1000    # maximum number of iterations
                                         , g_tol = 1e-8         # gradient tolerance
                                         , store_trace = true   # store the optimization trace
                                         , show_trace = true    # print progress to stdout
                                         , show_warnings = true  # show warnings
                                        )
                        ) 
    
    #TEST
    #params
    #init_gamma
    #res_gamma.minimizer
    #d0 = Gamma( res_gamma.minimizer[1],res_gamma.minimizer[2]) #(shape, scale)
    #d = truncated( d0, lower, upper)
    #test_d=rand(d,10000)
    #histogram(test_d,alpha=0.5, normalize=:pdf)
    #test_d0=rand(d0,10000)
    #histogram!(test_d0,alpha=0.5, normalize=:pdf, color=:green)
    #histogram!(times, normalize=:pdf, color=:red, alpha=0.5)
    # Overlay fitted (truncated) gamma PDF
    #xs = range(0, upper; length=300)
    #norm = cdf(d, upper) # normalization for truncation
    #tinf_pdf_vals = pdf.(d, xs) ./ norm # truncated PDF
    # Plot tinf fit
    #plot!(xs, tinf_pdf_vals; label="tinf truncated gamma fit", color=:blue, lw=2)
    #plot!(xs, pdf.(Gamma(18,3),xs); label="tinf truncated gamma fit", color=:blue, lw=2)

    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converged( res_gamma )
        nll_g = nll_trunc_gamma( Optim.minimizer( res_gamma ) )
        k = length( Optim.minimizer( res_gamma ) ) # Number of parameters estimated
        aic_g = 2*k + 2*nll_g 
        fit_results["Gamma"] = ( nll_g, aic_g, Optim.minimizer( res_gamma ) )
    #else
    #    fit_results["Gamma"] = ("NA", "NA", "NA")
    #end

    # Log-normal
    res_lognorm = optimize(nll_trunc_lognorm, init_lognorm, NelderMead()
                            , Optim.Options(
                                              iterations = 2000    # maximum number of iterations
                                            , g_tol = 1e-8         # gradient tolerance
                                            , store_trace = true   # store the optimization trace
                                            , show_trace = true    # print progress to stdout
                                            , show_warnings = true  # show warnings
                                            )
                            )
    
    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converged( res_lognorm ) #res_lognorm.exitflag == :Success
        nll_ln = nll_trunc_lognorm( Optim.minimizer( res_lognorm ) )
        k = length( Optim.minimizer( res_lognorm ) ) # Number of parameters estimated
        aic_ln = 2*k + 2*nll_ln
        fit_results["LogNormal"] = ( nll_ln, aic_ln, Optim.minimizer( res_lognorm ) )
    #else
    #    fit_results["LogNormal"] = ("NA", "NA", "NA")
    #end

    # Weibull
    res_weibull = optimize(nll_trunc_weibull, init_weibull, NelderMead()
                            , Optim.Options(
                                              iterations = 2000    # maximum number of iterations
                                            , g_tol = 1e-8         # gradient tolerance
                                            , store_trace = true   # store the optimization trace
                                            , show_trace = true    # print progress to stdout
                                            , show_warnings = true  # show warnings
                                            )
    )
    
    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converge( res_weibull ) #res_weibull.exitflag == :Success
        nll_w = nll_trunc_weibull(Optim.minimizer(res_weibull))
        k = length(Optim.minimizer(res_weibull)) # Number of parameters estimated
        aic_w = 2*k + 2*nll_w
        fit_results["Weibull"] = (nll_w, aic_w, Optim.minimizer(res_weibull))
    #else
    #   fit_results["Weibull"] = ("NA", "NA", "NA")
    #end

    # Negative Binomial
    res_nbinom = optimize(nll_trunc_nbinom, init_nbinom, NelderMead()
                            , Optim.Options(
                                              iterations = 2000    # maximum number of iterations
                                            , g_tol = 1e-8         # gradient tolerance
                                            , store_trace = true   # store the optimization trace
                                            , show_trace = true    # print progress to stdout
                                            , show_warnings = true  # show warnings
                                            )
    )

    # Check if optimisation of negative log-likelihood has converged
    # before savings results
    #if Optim.converge( res_nbinom ) #res_weibull.exitflag == :Success
        nll_nb = nll_trunc_nbinom( Optim.minimizer(res_nbinom) )
        k = length(Optim.minimizer(res_nbinom)) # Number of parameters estimated
        aic_nb = 2*k + 2*nll_nb
        fit_results["NegativeBinomial"] = (nll_nb, aic_nb, Optim.minimizer(res_nbinom))
    #else
    #    fit_results["NegativeBinomial"] = ("NA", "NA", "NA")
    #end


    # Report results
    #if isempty(fit_results)
    #    println("All fits failed.")
    #else
    #    println("Goodness-of-fit (AIC) for each distribution:")
    #    for (dist, (_nll, aic, params)) in results
    #        println("  $dist: AIC = $aic, parameters = $params")
    #    end
        # Select best distribution
    #    best_dist = findmin([(aic, dist) for (dist, (_nll, aic, _params)) in results])[2]
    #    println("Best-fitting distribution: $best_dist")
    #end
    return( [fit_results] )
end


#Testing gamma parameters so set initial params close to final params
#data = rand(truncated(Gamma(12.0, 100.0), 0, 60), 1000)
#data = rand(truncated(Weibull(20.0, 80.0), 0, 60), 1000)
#histogram(data)

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
#r = "TLI7"
#r = icu_regions[2]

    # Times to fit truncated gamma distribution to
    times = tinf_by_region_dict[ r ]
    #Test
    #histogram(times)

    # Define truncation range
    lower = 0 #0.00001
    upper = maxtime

    # Initial guesses
    #init_gamma = [10.0, 5.0] # shape*scale = mean. Median is ~50 for the tinf data # [1.0, 1.0]
    init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    init_weibull = [20.0, 80.0]

    fit_tinf_results = fit_multi_dist(; times
                                    , init_gamma, init_nbinom, init_lognorm, init_weibull
                                    , lower, upper)
    println(fit_tinf_results)

    # Fit truncated gamma distribution to tinf from earliest treport
    #init_params = [12.0, 6.0] #[1.0, 1.0]
    #res_gamma = optimize(negloglik_gamma, init_params, BFGS())
    #res_gamma = optimize(nll_trunc_gamma, init_params, BFGS())
    #shape_hat, scale_hat = Optim.minimizer(res_gamma)

    shape_hat, scale_hat = fit_tinf_results[1]["Gamma"][3]
    tinf_fitted_gamma = Gamma(shape_hat, scale_hat)
    tinf_trunc_gamma = truncated(tinf_fitted_gamma, 0.0, upper)

    # Plot histogram (normalized to probability density)
    p1 = histogram(times; bins=30
             , normalize=:pdf
             , label="tinf"
             , color=:lightblue, alpha=0.5, legend=:topleft
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
    tinf_median = quantile(times, 0.5) # TODO Includes zero values which were removed for stat dist fitting
    #tinf_median = median(times)
    p1 = vline!([tinf_median], color=:blue
            , label="tinf median")

    # Add median for fitted tinf
    tinf_median_fit = quantile(tinf_trunc_gamma, 0.5)
    p1 = vline!([tinf_median_fit], color=:blue, linestyle=:dash
            , label="tinf median truncated gamma fit")

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
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    init_weibull = [20.0, 80.0]

    fit_tinf_treport_results = fit_multi_dist(; times
                                               , init_gamma, init_nbinom, init_lognorm, init_weibull
                                               , lower, upper)
    println(fit_tinf_treport_results)

    #result = optimize(negloglik_gamma, init_params, BFGS())
    
    #shape_hat, scale_hat = Optim.minimizer(result)
    shape_hat, scale_hat = fit_tinf_treport_results[1]["Gamma"][3]
    tinf_treport_fitted_gamma = Gamma(shape_hat, scale_hat)
    tinf_treport_median_fit = median( tinf_treport_fitted_gamma )

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
median_times_by_region_sorted = sort( median_times_by_region , :treport_median_fit )
println( median_times_by_region_sorted )
CSV.write("scripts/median_TD_by_region_plots/TD_by_region_sorted.csv", median_times_by_region_sorted)

# TODO possibly resimulate
### Apply fitted statistical distribution to the minimum report time (treport / TD)

# Resimulate tsample and treport based on truncated fit to tinf

# Draw 1000 tinf from fitted distribution
tinf_fit_samples = rand(tinf_trunc_gamma, 1000)
histogram!(tinf_fit_samples, alpha = 0.1)
histogram!(tinf_fit_samples; bins=30, normalize=true, label="tinf fit samples", color=:lightgreen, alpha=0.6, legend=:topleft)

# Re-simulate tsample and treport


