### Similar to 'median_TD_by_region.jl' but looking at the variability 
### of the output with the number of simulations
### The input is an extract from the simulation output: G df filtered for ICU cases only

### Compute median time to detection (TD) for all ITL2 regions individually

#=
Outline
 1) Extract event timing information by region from simulation 
 2) Fit statistical distributions to event timing information 
 3) Report estimated TD by region
 4) Loop n times and collect estimated TD by region
 5) Repeat for different numbers of simulations
 6) Plot estimated TD by region range against number of simulations
 
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

using NBPMscape
using JLD2
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
using Pkg
Pkg.add("Optim")
using Optim
using RData 
using CSV 
using PowerAnalyses
using StatsPlots
#### Functions to fit to TD simulated distributions

# Function to compute the negative log-likelihood of a gamma distribution.
# This can then be optimised to fit truncated data (the tinf times which are truncated at the simulation maxtime).

# Negative log-likelihood for truncated Gamma
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

## Discretized version of gamma distribution
# Discretize the gamma distribution into a probability mass function (PMF) over specified bins
function discretize_gamma_pmf(shape::Float64, scale::Float64, bins::Vector{Float64})
    d = Gamma(shape, scale)
    cdfs = cdf.(d, bins)
    return diff(cdfs)
end

# Assign values to bins using searchsortedlast (faster than manual loop)
function assign_bins(values::Vector{Float64}, bins::Vector{Float64})
    return searchsortedlast.(Ref(bins), values) .- 1
end

# Negative log-likelihood for discretized gamma model
#function nll_disc_trunc_gamma(params::Vector{Float64}, times::Vector{Float64}, trunc::Tuple{Float64, Float64}; nbins::Int=100)
#function nll_disc_trunc_gamma(;params::Vector{Float64}, times::Vector{Float64}, trunc::Vector{Int, Int}, nbins::Int=100)
function nll_disc_trunc_gamma(;params, times, trunc, nbins)
    shape, scale = params
    if shape <= 0 || scale <= 0
        return 1e10
    end

    lower, upper = trunc
    bins = range(lower, stop=upper, length=nbins+1) |> collect

    pmf = discretize_gamma_pmf(shape, scale, bins)
    if any(pmf .<= 0)
        return 1e10
    end

    bin_indices = assign_bins(times, bins)
    bin_indices = clamp.(bin_indices, 1, nbins)  # Ensure indices are within bounds

    counts = countmap(bin_indices)
    ll = sum(counts[i] * log(pmf[i]) for i in keys(counts))

    return -ll
end

            # Negative log-likelihood for truncated Weibull
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

            
            # Function to compare distribution fits
            function fit_multi_dist(; times, init_gamma, init_disc_trunc_gamma, init_nbinom, init_lognorm, init_weibull, lower, upper)
                
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
                
                # Discretized gamma
                #params = init_disc_trunc_gamma, times = times, trunc = [0, 60], nbins = 60
                
                obj_f_disc_trunc_gamma = init_disc_trunc_gamma -> nll_disc_trunc_gamma(;params = init_disc_trunc_gamma, times = times, trunc = [0, 60], nbins = 60)
                #init_gamma=[16.466863607650694,2.8170337707969924]
                #res_gamma = optimize(nll_trunc_gamma, init_gamma, NelderMead() #GoldenSection()) #Brent()) #AcceleratedGradientDescent()) #MomentumGradientDescent()) #GradientDescent()) #ConjugateGradient()) #LBFGS()) #SimulatedAnnealing()) #NelderMead())
                res_disc_trunc_gamma = optimize(obj_f_disc_trunc_gamma, init_disc_trunc_gamma, NelderMead() #GoldenSection()) #Brent()) #AcceleratedGradientDescent()) #MomentumGradientDescent()) #GradientDescent()) #ConjugateGradient()) #LBFGS()) #SimulatedAnnealing()) #NelderMead())
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
                    nll_dg = nll_disc_trunc_gamma( params = Optim.minimizer( res_disc_trunc_gamma ), times = times, trunc = [lower, upper] , nbins = 60)
                    k = length( Optim.minimizer( res_disc_trunc_gamma ) ) # Number of parameters estimated
                    aic_dg = 2*k + 2*nll_dg 
                    fit_results["discGamma"] = ( nll_dg, aic_dg, Optim.minimizer( res_disc_trunc_gamma ) )
                #else
                #    fit_results["discGamma"] = ("NA", "NA", "NA")
                #end


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
                #if isempty(fit_results)
                #    println("All fits failed.")
                #else
                #    println("Goodness-of-fit (AIC) for each distribution:")
                #    for (dist, (_nll, aic, params)) in fit_results
                #        println("  $dist: negative log-likelihood = $_nll, AIC = $aic, parameters = $params")
                #    end
                # Select best distribution
                #best_dist = findmin([(aic, dist) for (dist, (_nll, aic, _params)) in fit_results])[1][2]
                #println("Best-fitting distribution: $best_dist ") #$fit_results[$("$best_dist")]")
                #end

                return( [fit_results] )

            end

# Load simulation
sims_G_icu_filter = load("covidlike-1.0-sims-regionentry_filtered_G_icu_combined.jld2", "sims_G_icu_filter")

# ITL2 region list
REGKEY.code

# How many replicates are in this simulation
nrep = length(sims_G_icu_filter)
# ICU sample proptions to test
psampled_test = [0.1, 0.25, 1.00]
# Number of times to repeat median TD measurement for each simulation
# in order to test variability
n_variabiality_test_repeats = 1 #10

# Repeat with different sampling proportion within region
for psampled in psampled_test

    # Create matrix to store the results
    # number of sim reps (rows) x regions (cols)
    # Median TD range by region
    n_rep_test_vals = 40000:5000:nrep # First batch = 2000:5000:nrep, where nrep = 37000
    nsimrep_range_fitted_median_td_range_by_region_variability = Array{Float64,2}(undef, length(n_rep_test_vals), length(REGKEY.code)) #
    
    # Repeat with different numbers of replicates
    # Start counter for number of repeat analyses
    n_simrep_test = 1

    for n_replicates in n_rep_test_vals #2000:5000:nrep
        #n_replicates = 60000
        #p_sampled = 0.25
        sims = sims_G_icu_filter[1:n_replicates]

        println("Number of replicates: $(n_replicates), ICU sample proportion: $(Int(psampled*100))%")

        # Create 3D arrays to store the results
        # sim reps (rows) x regions (cols) x repeats (slices)
        td_by_simrep_by_region_variability = Array{Float64,3}(undef, n_replicates, length(REGKEY.code), n_variabiality_test_repeats) #
        # median TD and median TD fit (rows) x regions (cols) x repeats (slices)
        td_by_simrep_by_region_variability_median = Array{Float64,3}(undef, 2, length(REGKEY.code), n_variabiality_test_repeats)
        
        # Repeat to see how result (median TD by region) varies
        for nvreps in 1:n_variabiality_test_repeats
            #nvreps = 1
            
            # Time to process sample and report results / declare detection
            turnaroundtime = 3

            # Create dictionary to hold the simulated time of infection for each ITL2 region
            tinf_by_region_dict         = Dict(name => Float64[] for name in REGKEY.code)
            tsample_by_region_dict      = Dict(name => Float64[] for name in REGKEY.code)
            tinf_tsample_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
            tinf_treport_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
            treport_by_region_dict      = Dict(name => Float64[] for name in REGKEY.code)
            # Create df to hold the number of ICU cases by region in each simulation replicate
            n_icu_by_region_by_simrep = DataFrame( zeros( length(sims) , length(REGKEY.code)), REGKEY.code )
            # Create df to record the TD by simulation replicate and by region
            td_by_simrep_by_region = DataFrame( [fill(Inf, length(sims)) for _ in 1:length(REGKEY.code)], REGKEY.code)
            
            # Loop through each simulation replicate and compile information
            for s in 1:length(sims)
                G = sims[s]
                if size(G,1) == 0
                    continue
                else
                    # Record number of ICU cases by region and by simulation
                    icu_count_by_region = countmap(G.homeregion)
                    
                    for (k, v) in icu_count_by_region
                        n_icu_by_region_by_simrep[ s, k ] = v
                    end
                    
                    # List regions with an ICU admission (i.e. possibility of getting a sample and detection)
                    icu_regions = unique(G.homeregion)
                    
                    # Loop through each region and record time of infection stats
                    for r in icu_regions #REGKEY.code
                        #r = icu_regions[1]
                        g_region = G[ (G.homeregion).==r, :]
                        
                        # Sample size 
                        n = rand( Binomial( size(g_region,1), psampled) ) 
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
                            tsample = map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region_sub) )
                            g_region_sub.tsample = tsample 

                            # Simulate reports times
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
                            
                            # Record TD (= min(treport)) by simulation replicate and by region
                            td_by_simrep_by_region[ s , r ] = min_treport
                        end
                    end            
                end

            end

            td_by_simrep_by_region_variability[:, :, nvreps] .= td_by_simrep_by_region
            
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
            
            # Plot to compare the median TD values by region and by simulation replicate
            median_td_by_region = [getindex( td_median_by_region, 1, i) for i in 1:size( td_median_by_region, 2)] 
            #println("Median of the median TD by region = ", round(median(median_td_by_region), digits=1))
            #println("Median of the minimum of regional TDs by sim rep = ", round( median(skipmissing(td_min_by_simrep)), digits=1))
            #println("with 95% CI: ",round(quantile( skipmissing(td_min_by_simrep) ,0.025), digits = 1)
            #        ," to ", round(quantile( skipmissing(td_min_by_simrep) ,0.975), digits=1))

            # Add median TD by region to results array in row 1
            # (fitted median TD by region will be added in row 2 below)
            td_by_simrep_by_region_variability_median[1, :, nvreps] = median_td_by_region
                        
            ### For each region, fit gamma distribution to infection times (tinf) corresponding the
            ### minimum report time (treport / TD)
            
            # Simulation maxtime
            maxtime = 60.0

            # Initialise df to record regions and the median tinf and treport (TD) times
            median_times_by_region = DataFrame( 
                ITL2_code = String[]
                , ITL2_name = String[]
                , tinf_median = Float64[]
                , tinf_median_gamma_fit = Float64[]
                , tinf_median_disc_gamma_fit = Float64[]
                , treport_median = Float64[]
                , treport_median_gamma_fit = Float64[]
                , treport_median_disc_gamma_fit = Float64[]
            )
            # Loop through all ITL2 regions

            for r in REGKEY.code # TEST: r = "TLI7"; r = icu_regions[2] ; r=REGKEY.code[1]
                #r = REGKEY.code[1]
                #r = REGKEY.code[27]
                
                # Times to fit truncated gamma distribution to
                times = tinf_by_region_dict[ r ]
                
                # Define truncation range
                lower = 0 #0.00001
                upper = maxtime

                # Initial guesses
                init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
                init_disc_trunc_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
                init_weibull = init_gamma #[20.0, 80.0]
                init_nbinom = [20.0, 0.5]
                init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
                
                fit_tinf_results = fit_multi_dist(; times
                                                , init_gamma, init_disc_trunc_gamma, init_nbinom, init_lognorm, init_weibull
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
                
                # discretized Gamma
                disc_gamma_shape_fit, disc_gamma_scale_fit = fit_tinf_results[1]["discGamma"][3]
                tinf_fitted_disc_gamma = Gamma(disc_gamma_shape_fit, disc_gamma_scale_fit)
                tinf_fitted_disc_gamma_trunc = truncated(tinf_fitted_disc_gamma, 0.0, upper)
                # Overlay fitted (truncated) gamma PDF
                xs = range(0, upper; length=300)
                norm_disc_gamma = cdf(tinf_fitted_disc_gamma, upper) # normalization for truncation
                tinf_pdf_vals_disc_gamma = pdf.(tinf_fitted_disc_gamma, xs) ./ norm_disc_gamma # truncated PDF
                tinf_pdf_vals_disc_gamma_trunc = pdf.(tinf_fitted_disc_gamma_trunc, xs) #./ norm_gamma # truncated PDF
                
                # Weibull
                weibull_shape_fit, weibull_scale_fit = fit_tinf_results[1]["Weibull"][3]
                tinf_fitted_weibull = Gamma(weibull_shape_fit, weibull_scale_fit)
                tinf_fitted_weibull_trunc = truncated(tinf_fitted_weibull, 0.0, upper)
                # Overlay fitted (truncated) gamma PDF
                norm_weibull = cdf(tinf_fitted_weibull, upper) # normalization for truncation
                tinf_pdf_vals_weibull = pdf.(tinf_fitted_weibull, xs) ./ norm_weibull # truncated PDF
                tinf_pdf_vals_weibull_trunc = pdf.(tinf_fitted_weibull_trunc, xs) #./ norm_weibull # truncated PDF
                
                # Add median for tinf
                tinf_median = quantile(times, 0.5)
                times_no_zeros = sort( filter(x -> x != 0, times) )
                tinf_median_no_zeros = quantile(times_no_zeros, 0.5)
                
                # Add median for fitted tinf
                tinf_median_fit_gamma = quantile(tinf_fitted_gamma_trunc, 0.5)
                tinf_median_fit_disc_gamma = quantile(tinf_fitted_disc_gamma_trunc, 0.5)
                tinf_median_fit_weibull = quantile(tinf_fitted_weibull_trunc, 0.5)
                
                # Add median for treport
                treport_median = quantile(treport_by_region_dict[ r ], 0.5)
                               
                # Fit gamma distribution to (treport - tinf) from earliest treport
                times = tinf_treport_by_region_dict[ r ]
                
                # Initial guesses
                #init_params = [1.0, 1.0]
                init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
                init_disc_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
                init_weibull = init_gamma
                init_nbinom = [20.0, 0.5]
                init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
                
                fit_tinf_treport_results = fit_multi_dist(; times
                                                        , init_gamma, init_disc_trunc_gamma, init_nbinom, init_lognorm, init_weibull
                                                        , lower, upper)
                
                # Gamma fit
                tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit = fit_tinf_treport_results[1]["Gamma"][3]
                tinf_treport_fitted_gamma = Gamma(tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit)
                tinf_treport_median_gamma_fit = median( tinf_treport_fitted_gamma ) # Note that this distribution is not truncated
                
                # Overlay fitted gamma PDF
                tinf_treport_pdf_vals_gamma = pdf.(tinf_treport_fitted_gamma, xs)
                
                # discretized Gamma fit
                tinf_treport_disc_gamma_shape_fit, tinf_treport_disc_gamma_scale_fit = fit_tinf_treport_results[1]["discGamma"][3]
                tinf_treport_fitted_disc_gamma = Gamma(tinf_treport_disc_gamma_shape_fit, tinf_treport_disc_gamma_scale_fit)
                tinf_treport_median_disc_gamma_fit = median( tinf_treport_fitted_disc_gamma ) # Note that this distribution is not truncated
                
                # Overlay fitted discretized gamma PDF
                tinf_treport_pdf_vals_disc_gamma = pdf.(tinf_treport_fitted_disc_gamma, xs)
                
                # Add median for (treport - tinf)
                tinf_treport_median = quantile(tinf_treport_by_region_dict[ r ], 0.5) # Note does not include any zeros
                
                # Add median for fitted (treport - tinf)
                tinf_treport_median_fit_gamma = quantile(tinf_treport_fitted_gamma, 0.5)
                tinf_treport_median_fit_disc_gamma = quantile(tinf_treport_fitted_disc_gamma, 0.5)
                
                # Add row to df to record the median value by region
                new_row = ( ITL2_code = r
                            , ITL2_name = REGKEY[REGKEY.code .== r, :name][1]
                            , tinf_median = tinf_median
                            , tinf_median_gamma_fit = tinf_median_fit_gamma
                            , tinf_median_disc_gamma_fit = tinf_median_fit_disc_gamma
                            , treport_median = treport_median
                            , treport_median_gamma_fit = (tinf_treport_median_fit_gamma + tinf_median_fit_gamma)
                            , treport_median_disc_gamma_fit = (tinf_treport_median_fit_disc_gamma + tinf_median_fit_disc_gamma)
                            )
                push!( median_times_by_region, new_row )
                
            end

            #println( median_times_by_region )
            #median_times_by_region_sorted = sort( median_times_by_region , :treport_median_fit )
            #println( median_times_by_region_sorted )

                    
            # Add median and fit median to array
            td_by_simrep_by_region_variability_median[1,:,nvreps] = median_times_by_region.treport_median
            td_by_simrep_by_region_variability_median[2,:,nvreps] = median_times_by_region.treport_median_gamma_fit #treport_median_fit

        end # of loop for test of variability of median TD by region'/
        
        # Save 3D array containing sim reps (rows) x regions (cols) x repeats (slices)
        @save "scripts/median_TD_by_region/sim_regionentry/variability_check/td_by_simrep_by_region_variability_$(n_replicates)_$(Int(psampled*100)).jld2" td_by_simrep_by_region_variability
        @save "scripts/median_TD_by_region/sim_regionentry/variability_check/td_by_simrep_by_region_variability_median_$(n_replicates)_$(Int(psampled*100)).jld2" td_by_simrep_by_region_variability_median
        
        # Compute range in median TD (simulated and fitted)
        m, n, _ = size(td_by_simrep_by_region_variability_median)
        max_median_TD_by_region = [maximum(td_by_simrep_by_region_variability_median[i, j, :]) for i in 1:m, j in 1:n]
        min_median_TD_by_region = [minimum(td_by_simrep_by_region_variability_median[i, j, :]) for i in 1:m, j in 1:n]
        range_median_TD_by_region = max_median_TD_by_region - min_median_TD_by_region
        
        #plot(range_median_TD_by_region[1,:])
        #plot!(range_median_TD_by_region[2,:])
        
        # Record range of median fitted TD by region with a row for each number of simulations
        nsimrep_range_fitted_median_td_range_by_region_variability[n_simrep_test,:] = range_median_TD_by_region[2,:]
        
        n_simrep_test = n_simrep_test + 1

    end # of loop for different number of simulations

    # Save 2D array containing sim reps (rows) x regions (cols) x repeats (slices)
    @save "scripts/median_TD_by_region/sim_regionentry/variability_check/nsimrep_range_fitted_median_td_range_by_region_variability_$(Int(psampled*100)).jld2" nsimrep_range_fitted_median_td_range_by_region_variability
        
end # of loop different psampled values

### View variability results
### Load matrices containing n_sim_reps (rows) x regions (cols) = range of fitted median TDs
nsimrep_range_fitted_median_td_range_by_region_variability_100 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/nsimrep_range_fitted_median_td_range_by_region_variability_100.jld2", "nsimrep_range_fitted_median_td_range_by_region_variability")
nsimrep_range_fitted_median_td_range_by_region_variability_25 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/nsimrep_range_fitted_median_td_range_by_region_variability_25.jld2", "nsimrep_range_fitted_median_td_range_by_region_variability")
nsimrep_range_fitted_median_td_range_by_region_variability_10 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/nsimrep_range_fitted_median_td_range_by_region_variability_10.jld2", "nsimrep_range_fitted_median_td_range_by_region_variability")

plot(permutedims(nsimrep_range_fitted_median_td_range_by_region_variability_100[8:8,:]))
plot!(permutedims(nsimrep_range_fitted_median_td_range_by_region_variability_25[8:8,:]))
plot!(permutedims(nsimrep_range_fitted_median_td_range_by_region_variability_10[8:8,:]))

size(nsimrep_range_fitted_median_td_range_by_region_variability_100)

scatter([2,7,12,17,22,27,32,37],nsimrep_range_fitted_median_td_range_by_region_variability_100[1:8,:])
scatter([2,7,12,17,22,27,32,37],nsimrep_range_fitted_median_td_range_by_region_variability_25[1:8,:])
scatter([2,7,12,17,22,27,32,37],nsimrep_range_fitted_median_td_range_by_region_variability_10[1:8,:])

scatter([2,7,12,17,22,27,32,37],mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_100[1:8,:]; dims=2))
scatter!([2,7,12,17,22,27,32,37],mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_25[1:8,:];  dims=2))
scatter!([2,7,12,17,22,27,32,37],mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_10[1:8,:];  dims=2))

TD_range_median_p100 = mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_100[:,:]; dims=2)
TD_range_median_p25 = mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_25[:,:]; dims=2)
TD_range_median_p10 = mapslices(median, nsimrep_range_fitted_median_td_range_by_region_variability_10[:,:]; dims=2)
TD_range_min_p100 = mapslices(minimum, nsimrep_range_fitted_median_td_range_by_region_variability_100[:,:]; dims=2)
TD_range_min_p25 = mapslices(minimum, nsimrep_range_fitted_median_td_range_by_region_variability_25[:,:]; dims=2)
TD_range_min_p10 = mapslices(minimum, nsimrep_range_fitted_median_td_range_by_region_variability_10[:,:]; dims=2)
TD_range_max_p100 = mapslices(maximum, nsimrep_range_fitted_median_td_range_by_region_variability_100[:,:]; dims=2)
TD_range_max_p25 = mapslices(maximum, nsimrep_range_fitted_median_td_range_by_region_variability_25[:,:]; dims=2)
TD_range_max_p10 = mapslices(maximum, nsimrep_range_fitted_median_td_range_by_region_variability_10[:,:]; dims=2)

p1 = plot( [2,7,12,17,22,27,32,37]
    , TD_range_median_p10 
    , linewidth = 4
    #, ribbon = (TD_range_min_p10, TD_range_max_p10)
    , label = "psampled = 0.10"
    ,xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , ylimit = [0,2]
    )

p1 = plot!( [2,7,12,17,22,27,32,37]
    , TD_range_median_p25 
    , linewidth = 4
    #, ribbon = (TD_range_min_p25, TD_range_max_p25)
    , label = "psampled = 0.25"
    #, ylimit = [0,7]
    )

p1 = plot!( [2,7,12,17,22,27,32,37]
    , TD_range_median_p100
    , linewidth = 4
    #, ribbon = (TD_range_min_p100, TD_range_max_p100)
    , label = "psampled = 1.00"
    #, ylimit = [0,7]
    )

p2 = plot( [2,7,12,17,22,27,32,37]
    , TD_range_median_p10 
    , linewidth = 4
    , ribbon = (TD_range_min_p10, TD_range_max_p10)
    , label = "psampled = 0.10"
    , xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , ylimit = [0,7]
    )

p3 = plot( [2,7,12,17,22,27,32,37]
    , TD_range_median_p25 
    , linewidth = 4
    , ribbon = (TD_range_min_p25, TD_range_max_p25)
    , label = "psampled = 0.25"
    , xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , color = "red"
    , ylimit = [0,7]
    )

p4 = plot( [2,7,12,17,22,27,32,37]
    , TD_range_median_p100
    , linewidth = 4
    , ribbon = (TD_range_min_p100, TD_range_max_p100)
    , label = "psampled = 1.00"
    , xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , color = "green"
    , ylimit = [0,7]
    )

plot( size=(800, 1200), p1, p2, p3, p4, layout = (4,1))
savefig("scripts/median_TD_by_region/sim_regionentry/variability_check/median_TD_variability_10_repeats.png")


### Merge the two different variability analysis results sims_test
# This is the difference between the maximum and minimum of the fitted median TDs
# across the n repeat analyses, for each region
x_axis_sim_reps = [2,7,12,17,22,27,32,37,40,45,50,55,60] # These relate to the rows in the nsimrep_range_fitted_median_td_range_by_region_variability_$(Int(psampled*100)).jld2 files

# Read in first analysis of lower sim rep numbers (up to 37,000) and
# second analysis (up to 60,000)
# Low sim rep analysis
low_simrep_variability_p10 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_2000_37000/nsimrep_range_fitted_median_td_range_by_region_variability_10.jld2"
, "nsimrep_range_fitted_median_td_range_by_region_variability")                                

low_simrep_variability_p25 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_2000_37000/nsimrep_range_fitted_median_td_range_by_region_variability_25.jld2"
                                , "nsimrep_range_fitted_median_td_range_by_region_variability")
low_simrep_variability_p100 = load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_2000_37000/nsimrep_range_fitted_median_td_range_by_region_variability_100.jld2"
                                , "nsimrep_range_fitted_median_td_range_by_region_variability")
# High sim rep analysis
high_simrep_variability_p10 =load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_40000_60000/nsimrep_range_fitted_median_td_range_by_region_variability_10.jld2"
                                , "nsimrep_range_fitted_median_td_range_by_region_variability")
high_simrep_variability_p25 =load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_40000_60000/nsimrep_range_fitted_median_td_range_by_region_variability_25.jld2"
                                , "nsimrep_range_fitted_median_td_range_by_region_variability")
high_simrep_variability_p100 =load("scripts/median_TD_by_region/sim_regionentry/variability_check/summary_40000_60000/nsimrep_range_fitted_median_td_range_by_region_variability_100.jld2"
                                , "nsimrep_range_fitted_median_td_range_by_region_variability")
# Merge
diff_median_TD_10_repeats_p10 = vcat(low_simrep_variability_p10, high_simrep_variability_p10)
diff_median_TD_10_repeats_p25 = vcat(low_simrep_variability_p25, high_simrep_variability_p25)
diff_median_TD_10_repeats_p100 = vcat(low_simrep_variability_p100, high_simrep_variability_p100)

# Add column names and row names, while also transposing to rows = regions and cols = number of sim reps
# p10
diff_median_TD_10_repeats_p10_df = DataFrame( transpose(diff_median_TD_10_repeats_p10)
                                            , Symbol.(string.(x_axis_sim_reps)))
diff_median_TD_10_repeats_p10_df.ITL2_code = REGKEY.code
#diff_median_TD_10_repeats_p10_df.ITL2_name = REGKEY.code
CSV.write("scripts/median_TD_by_region/sim_regionentry/variability_check/diff_median_TD_10_repeats_p10.csv"
        , diff_median_TD_10_repeats_p10_df)
# p25
diff_median_TD_10_repeats_p25_df = DataFrame( transpose(diff_median_TD_10_repeats_p25)
                                            , Symbol.(string.(x_axis_sim_reps)))
diff_median_TD_10_repeats_p25_df.ITL2_code = REGKEY.code
#diff_median_TD_10_repeats_p25_df.ITL2_name = REGKEY.code
CSV.write("scripts/median_TD_by_region/sim_regionentry/variability_check/diff_median_TD_10_repeats_p25.csv"
        , diff_median_TD_10_repeats_p25_df)
# p100
diff_median_TD_10_repeats_p100_df = DataFrame( transpose(diff_median_TD_10_repeats_p100)
                                            , Symbol.(string.(x_axis_sim_reps)))
diff_median_TD_10_repeats_p100_df.ITL2_code = REGKEY.code
#diff_median_TD_10_repeats_p100_df.ITL2_name = REGKEY.code
CSV.write("scripts/median_TD_by_region/sim_regionentry/variability_check/diff_median_TD_10_repeats_p100.csv"
        , diff_median_TD_10_repeats_p100_df)


### Plot
TD_range_median_p100 = mapslices(median, diff_median_TD_10_repeats_p100[:,:]; dims=2)
TD_range_median_p25 = mapslices(median, diff_median_TD_10_repeats_p25[:,:]; dims=2)
TD_range_median_p10 = mapslices(median, diff_median_TD_10_repeats_p10[:,:]; dims=2)
TD_range_min_p100 = mapslices(minimum, diff_median_TD_10_repeats_p100[:,:]; dims=2)
TD_range_min_p25 = mapslices(minimum, diff_median_TD_10_repeats_p25[:,:]; dims=2)
TD_range_min_p10 = mapslices(minimum, diff_median_TD_10_repeats_p10[:,:]; dims=2)
TD_range_max_p100 = mapslices(maximum, diff_median_TD_10_repeats_p100[:,:]; dims=2)
TD_range_max_p25 = mapslices(maximum, diff_median_TD_10_repeats_p25[:,:]; dims=2)
TD_range_max_p10 = mapslices(maximum, diff_median_TD_10_repeats_p10[:,:]; dims=2)

p1 = plot( x_axis_sim_reps
    , TD_range_median_p10 
    , linewidth = 4
    #, ribbon = (TD_range_min_p10, TD_range_max_p10)
    , label = "psampled = 0.10"
    ,xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , ylimit = [0,2]
    )

p1 = plot!( x_axis_sim_reps
    , TD_range_median_p25 
    , linewidth = 4
    #, ribbon = (TD_range_min_p25, TD_range_max_p25)
    , label = "psampled = 0.25"
    #, ylimit = [0,7]
    )

p1 = plot!( x_axis_sim_reps
    , TD_range_median_p100
    , linewidth = 4
    #, ribbon = (TD_range_min_p100, TD_range_max_p100)
    , label = "psampled = 1.00"
    #, ylimit = [0,7]
    )

p2 = plot( x_axis_sim_reps
    , TD_range_median_p10 
    , linewidth = 4
    , ribbon = (TD_range_min_p10, TD_range_max_p10)
    , label = "psampled = 0.10"
    , xlabel = "sim replicates ('000s)"
    , ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , ylimit = [0,7]
    )

p3 = plot( x_axis_sim_reps
    , TD_range_median_p25 
    , linewidth = 4
    , ribbon = (TD_range_min_p25, TD_range_max_p25)
    , label = "psampled = 0.25"
    , xlabel = "sim replicates ('000s)"
    #, ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , color = "red"
    , ylimit = [0,7]
    )

p4 = plot( x_axis_sim_reps
    , TD_range_median_p100
    , linewidth = 4
    , ribbon = (TD_range_min_p100, TD_range_max_p100)
    , label = "psampled = 1.00"
    , xlabel = "sim replicates ('000s)"
    #, ylabel = "Max-Min median TD\n for individual ITL2 regions\n over 10 repeat analyses"
    , color = "green"
    , ylimit = [0,7]
    )

plot( size=(1200, 800), p1, p2, p3, p4, layout = (1,4))
savefig("scripts/median_TD_by_region/sim_regionentry/variability_check/median_TD_variability_10_repeats_nrep60k.png")

# Also violin plots because the median is very low and close to the minimum so
# probably very few regions at the upper end

function violin_plot(;sim_reps, median_diff_data, plot_col, p_samp_val, y_max)
    for k in 1:size(median_diff_data,1)
        if k == 1
            plot_output = violin( fill(sim_reps[k],length(median_diff_data[k,:]))
                                , median_diff_data[k,:]
                                , xlabel = "sim replicates ('000s)"
                                , ylabel = "Max-Min TD\n for individual ITL2 regions\n over 10 repeat analyses"
                                , title = "psampled = $(p_samp_val)" #0.10
                                , legend = false
                                , ylimit = [0,y_max]#[0,7]
                                , color = plot_col #"blue"
                                )
        else
        plot_output = violin!( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p10[k,:]))
                                , diff_median_TD_10_repeats_p10[k,:]  )
        end
    end
    return(plot_output)
end
maximum(diff_median_TD_10_repeats_p10[1,:])

pp1 = violin_plot(sim_reps = x_axis_sim_reps
                , median_diff_data = diff_median_TD_10_repeats_p10
                , plot_col = "blue"
                , p_samp_val = 0.10
                , y_max = 7)
display(pp1)

# psampled = 10
for k in 1:size(diff_median_TD_10_repeats_p10,1)
    if k == 1
        pp1 = violin( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p10[k,:]))
                    , diff_median_TD_10_repeats_p10[1,:]
                    , xlabel = "sim replicates ('000s)"
                    , ylabel = "Max-Min TD\n for individual ITL2 regions\n over 10 repeat analyses"
                    , title = "psampled = 0.10"
                    , legend = false
                    , ylimit = [0,7]
                    , color = "blue"
                    )
    else
    pp1 = violin!( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p10[k,:]))
                , diff_median_TD_10_repeats_p10[k,:]
                , color = "blue"  )
    end
end
display(pp1)

# psampled = 25
for k in 1:size(diff_median_TD_10_repeats_p25,1)
    if k == 1
        pp2 = violin( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p25[k,:]))
                    , diff_median_TD_10_repeats_p25[1,:]
                    , xlabel = "sim replicates ('000s)"
                    , ylabel = "Max-Min TD\n for individual ITL2 regions\n over 10 repeat analyses"
                    , title = "psampled = 0.25"
                    , legend = false
                    , ylimit = [0,7]
                    , color = "red"
                    )
    else
    pp2 = violin!( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p25[k,:]))
                , diff_median_TD_10_repeats_p25[k,:]
                , color = "red"  )
    end
end
display(pp2)
maximum(diff_median_TD_10_repeats_p10[1,:])

# psampled = 100
for k in 1:size(diff_median_TD_10_repeats_p100,1)
    if k == 1
        pp3 = violin( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p100[k,:]))
                    , diff_median_TD_10_repeats_p100[1,:]
                    , xlabel = "sim replicates ('000s)"
                    , ylabel = "Max-Min TD\n for individual ITL2 regions\n over 10 repeat analyses"
                    , title = "psampled = 1.00"
                    , legend = false
                    , ylimit = [0,7]
                    , color = "green"
                    )
    else
    pp3 = violin!( fill(x_axis_sim_reps[k],length(diff_median_TD_10_repeats_p100[k,:]))
                , diff_median_TD_10_repeats_p100[k,:]
                , color = "green"  )
    end
end
display(pp3)

plot( size=(1200, 800), pp1, pp2, pp3, layout = (1,3))
savefig("scripts/median_TD_by_region/sim_regionentry/variability_check/median_TD_variability_10_repeats_nrep60k_violin.png")

### Plot with maximum number of replicates
ppp1 = violin( fill( 0.1 , length(diff_median_TD_10_repeats_p10[13,:]))
                    , diff_median_TD_10_repeats_p10[13,:]
                    , xlabel = "ICU sampling proportion"
                    , ylabel = "Max-Min TD\n for individual ITL2 regions\n over 10 repeat analyses"
                    , title = "Number of sim replicates = 60,000"
                    , label = "psampled = 0.10"
                    , legend = true
                    , ylimit = [0,1]
                    , color = "blue"
                    , alpha = 0.5
                    )
ppp1 = violin!( fill( 0.25 , length(diff_median_TD_10_repeats_p25[13,:]))
                    , diff_median_TD_10_repeats_p25[13,:]
                    , label = "psampled = 0.25"
                    , ylimit = [0,1]
                    , color = "red"
                    , alpha = 0.5)

ppp1 = violin!( fill( 1.00 , length(diff_median_TD_10_repeats_p100[13,:]))
                    , diff_median_TD_10_repeats_p100[13,:]
                    , label = "psampled = 1.00"
                    , ylimit = [0,1]
                    , color = "green"
                    , alpha = 0.5
                    )

savefig("scripts/median_TD_by_region/sim_regionentry/variability_check/median_TD_variability_10_repeats_nrep60k_only_violin.png")

### Variation in correlation

# 1 Start with analysis using max number of sim reps (=60,000)
# 2 Compute correlation matrix for fitted median TD for each region, 
#   with a separate correlation matrix for each repeat analysis (so 10 corr matrices)
# 3 Display three heatmaps:
#   3.1 median correlation (for each region across 10 repeat analyses)
#   3.2 absolute difference between min and max correlation (for each region across 10 repeat analyses)
#   3.3 % difference between min and max (for each region across 10 repeat analyses)

# Load file containing data for 60,000 replicates (the max so far) and 
td_by_simrep_by_region_variability = load("scripts/median_TD_by_region/sim_regionentry/variability_check/td_by_simrep_by_region_variability_60000_25.jld2", "td_by_simrep_by_region_variability")


####### BELOW UNUSED ################
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
            
            corr_df = DataFrame(corrmat, Symbol.(names(td_by_simrep_by_region)))
            maximum(skipmissing(corrmat))
            maximum(x for x in skipmissing(corrmat) if x != 1.0)
            minimum(skipmissing(corrmat))
            
            # Save correlation matrix to a file
            corr_df_w_names = insertcols!(corr_df, 1, :ITL2_region => y_names )
            #println(corr_df_w_names)
            #CSV.write("scripts/median_TD_by_region/sim_regionentry/$(n_replicates)_replicates/psampled_$(Int(psampled*100))/regional_TD_correlation_matrix_psampled_$(Int(psampled*100)).csv", corr_df_w_names)

            ## Plot a correlation matrix with a reduced number of ITL2 regions
            region_names = names(corr_df)
            subset_regions = ["TLI4","TLI3","TLI6","TLH4","TLG3","TLF2","TLF1","TLC2","TLE3","TLE4","TLD3","TLD7","TLJ3","TLJ1","TLK1","TLL2","TLM1","TLN0"]
            subset_regions = sort(subset_regions)
            #["TLC2", "TLD3", "TLD7", "TLE3", "TLE4", "TLF1", "TLF2", "TLG3", "TLH4"
            #, "TLI3", "TLI4", "TLI6", "TLJ1", "TLJ3", "TLK1", "TLL2", "TLM1", "TLN0"]
            # Subset regions not currently in region names - due to change of code/region
            missing_regions = subset_regions[ findall(i -> isnothing(i), indexin(subset_regions, region_names)) ]
            # TLC2 changed in 2025 - replace with TLC4
            # TLK1 changed in 2025 - replace with TLK5
            # TLL2 changed in 2025 - replace with TLL5
            subset_regions = [["TLC4","TLK5","TLL5"];subset_regions]
            subset_regions = sort(subset_regions)
            # Find indices of the subset regions in full list of region names
            subset_idx = findall(x -> x in subset_regions, region_names)
            # Subset the correlation matrix
            filtered_corr = corr_df_w_names[subset_idx.-1, [1;subset_idx]]
            
            # Save heatmap plot and data
            CSV.write("scripts/median_TD_by_region/sim_regionentry/$(n_replicates)_replicates/psampled_$(Int(psampled*100))/regional_TD_correlation_matrix_psampled_$(Int(psampled*100))_reduced.csv", filtered_corr)    




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