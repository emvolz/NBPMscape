#=

Summary
        Set of functions that can be used to fit statistical distributions to truncated 
        (e.g. by the max simulation time (maxtime) used in simtree and simforest) and/or 
        sparse (e.g. regional rather than national sampling) time to detection (TD) datasets. 
        The fitted statistical distribution can then be used to estimate the median time to detection.

List of functions in this file:

        - nll_trunc_gamma       Computes negative log-likelihood for truncated Gamma distribution fit
        - nll_trunc_normal      Computes negative log-likelihood for truncated Normal distribution fit
        - nll_trunc_weibull     Computes negative log-likelihood for truncated Weibull distribution fit
        - discretize_gamma_pmf  Discretize the Gamma distribution into a probability mass function (PMF) over 
                                specified bins.
        - assign_bins           Assign values to bins using searchsortedlast (faster than manual loop)
                                Used in function to fit a discretized Gamma distribution (nll_disc_gamma) to
                                assign continuous time values to discrete bins
        - nll_disc_gamma        Negative log-likelihood for discretized gamma model

        - dist_fit_plot     Plots data and fitted distribution
        - fit_multi_dist    Function to compare results for fitting Gamma and Weibull statistical distributions to data
                            Note that this could potentially be extended to include other statistical distributions.
        - erlang_truncated_means    Splits an Erlang distribution into two parts and returns the mean x value of each part

Description of regional time to detection (TD) estimation and requirement for distribution fitting:
        Extract region specific information from each forest in the outbreak simulation. 
        However, not all forests will contain ICU infections (able to be sampled) in a 
        particular ITL2 region. In that case, a robust way of estimating it would be to fit a
        parametric distribution (e.g. Gamma) to the TD's. 
        The simulation has maximum run time (maxtime) and so truncation needs to be accounted for.
        Accounting for truncation is more difficult for treport and so better to use tinf as 
        this is the parameter that is actually truncated by maxtime.
        Possible to fit to tinf and then fit separately to (treport - tinf), which is not truncated.
        These can then be summed to obtain an estimate of treport (= tinf + (treport - tinf))
        Finally report median estimates for TD (earliest treport) by region.

=#


"""
Function:       nll_trunc_gamma

Description:    Function to compute the negative log-likelihood of a gamma distribution.
                This can then be optimised to fit truncated data (e.g. the tinf times which are truncated at the simulation maxtime).

Arguments:      params  Two element vector containing the intial shape and scale parameters for the Gamma distribution
                times   Vector of time values to be fit with statistical distribution
                trunc   Two element vector containing the lower and upper limits of the truncation, e.g. [0, 90]

Returns:        Negative log-likelihood of the statistical distribution parameters (params) fitting the data (times)

Examples:   
                times = [0.09, 0.98, 1.21, 3.82, 2.77, 1.53, 1.57, 1.35, 1.51, 0.77]
                init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)]
                nll_trunc_gamma(;params = init_gamma, times = times, trunc = [0, 4])

"""
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
          
    return isfinite(-ll) ? -ll : 1e10
end


"""
Function:   discretize_gamma_pmf

Description:    Discretized version of Gamma distribution. 
                Discretize the Gamma distribution into a probability mass function (PMF) over 
                specified bins.

Arguments:      shape::Float64          Parameter for continuous Gamma distribution
                scale::Float64          Parameter for continuous Gamma distribution
                bins::Vector{Float64}   Number of bins in discretisation

Returns:        Vector representing values of defined continuous Gamma distribution discretized into the input bin edges.
                See example below.

Examples:       dgp = discretize_gamma_pmf(shape = 344.0, scale = 0.26
                                          , bins = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]
                                          )
                # Returns
                20-element Vector{Float64}:
                [ 3.323997277463239e-292, 5.630823812300795e-197, 1.0052059782315457e-144, 4.564559290809724e-110, 4.7519453236107885e-85, 3.964650110552051e-66
                 , 2.0600437488870638e-51, 8.966581051882102e-40, 1.7493619235719315e-30, 4.813530073250704e-23, 4.2372779372367254e-17, 2.189121695892049e-12
                 , 1.0553341009344165e-8, 6.834250680461053e-6, 0.0007970911273459353, 0.0212590215667613, 0.1571497039548055, 0.3740268022979566, 0.3212354494302465
                 , 0.10881266983516491]
                
"""
function discretize_gamma_pmf(; shape::Float64, scale::Float64, bins::Vector{Float64})
    d = Gamma(shape, scale)
    cdfs = cdf.(d, bins)
    return diff(cdfs)
end


"""
Function:   assign_bins

Description:    Assign values to bins using searchsortedlast (faster than manual loop)
                Used in function to fit a discretized Gamma distribution (nll_disc_gamma) to
                assign continuous time values to discrete bins

Arguments:  values::Vector{Float64}
            bins::Vector{Float64}

Returns:

Examples:
           assign_bins( 
                        vals = [ 3.323997277463239e-292, 5.630823812300795e-197, 1.0052059782315457e-144, 4.564559290809724e-110, 4.7519453236107885e-85, 3.964650110552051e-66
                                 , 2.0600437488870638e-51, 8.966581051882102e-40, 1.7493619235719315e-30, 4.813530073250704e-23, 4.2372779372367254e-17, 2.189121695892049e-12
                                 , 1.0553341009344165e-8, 6.834250680461053e-6, 0.0007970911273459353, 0.0212590215667613, 0.1571497039548055, 0.3740268022979566, 0.3212354494302465
                                 , 0.10881266983516491]
                        , bins = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]
                        )

"""
#function assign_bins(; vals::Vector{Float64}, bins::Vector{Float64})
function assign_bins(; vals, bins)
    return searchsortedlast.(Ref(bins), vals) .- 1
end

 

"""
Function:   nll_disc_gamma

Description:    Negative log-likelihood for discretized gamma model

Arguments:      params::Vector{Float64} Parameters for continuous Gamma distribution [shape, scale]
                times::Vector{Float64}  Data to fit discretized Gamma distribution to
                trunc::Vector           [lower, upper] limits of trucated distribution
                nbins::Int              Number of bins for discrisation

Returns:    Negative log-likelihood of discretized Gamma distribution fit to the data (times)

Examples:
            nll_disc_gamma(  params = [344.0, 0.26] # Gamma[shape, scale]
                           , times = [83.40972030623539, 85.16094984699521, 87.24240138225899, 84.83907644798674, 82.00932482461725, 83.483348202937, 88.66115627528256, 87.57949422080628, 89.61445158171773, 85.57725725192338]
                           , trunc = [0,100]
                           , nbins = 20
                           )
            # Returns 36.49

"""
function nll_disc_gamma(; params::Vector{Float64}, times::Vector{Float64}, trunc::Vector, nbins::Int=100) #trunc::Tuple{Float64, Float64}, nbins::Int=100)
    
    # Define Gamma distribution parameters
    shape, scale = params
    # Validate parameters
    if shape <= 0 || scale <= 0
        return 1e10
    end

    # Define range limits
    lower, upper = trunc
    # Define bin edges
    bins = range(lower, stop=upper, length=nbins+1) |> collect

    # Compute Gamma probability mass function values
    pmf = NBPMscape.discretize_gamma_pmf(shape = shape, scale = scale, bins = bins)
    # Validate values
    if any(pmf .<= 0)
        return 1e10
    end

    # Assign times to bins
    bin_indices = assign_bins(vals = times, bins = bins)
    bin_indices = clamp.(bin_indices, 1, nbins)  # Ensure indices are within bounds

    # Count the number of times in each bin
    counts = StatsBase.countmap(bin_indices)
    # And compute log-likelihood
    ll = sum(counts[i] * log(pmf[i]) for i in keys(counts))

    # Return negative log-likelihood
    return isfinite(-ll) ? -ll : 1e10
end


"""
Function:       nll_trunc_weibull

Description:    Function to compute the negative log-likelihood of a Weibull distribution fit.
                This can then be optimised to fit truncated data (e.g. the tinf times which are truncated at the simulation maxtime).

Arguments:      params  Two element vector containing the intial alpha and theta parameters for the Weibull distribution
                times   Vector of time values to be fit with statistical distribution
                trunc   Two element vector containing the lower and upper limits of the truncation, e.g. [0, 90]

Returns:        Negative log-likelihood of the statistical distribution parameters (params) fitting the data (times)

Examples:   
                times = [0.09, 0.98, 1.21, 3.82, 2.77, 1.53, 1.57, 1.35, 1.51, 0.77]
                init_weibull = [mean(times)^2 / var(times), var(times)/mean(times)]
                nll_trunc_weibull(;params = init_weibull, times = times, trunc = [0, 4])
"""
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
    
    return isfinite(-ll) ? -ll : 1e10
end


"""
Function:       nll_trunc_normal

Description:    Function to compute the negative log-likelihood of a Normal distribution fit.
                This can then be optimised to fit truncated data (e.g. the tinf times which are truncated at the simulation maxtime).

Arguments:      params  Two element vector containing the intial mu and sigma parameters for the Normal distribution
                times   Vector of time values to be fit with statistical distribution
                trunc   Two element vector containing the lower and upper limits of the truncation, e.g. [0, 90]

Returns:        Negative log-likelihood of the statistical distribution parameters (params) fitting the data (times)

Examples:   
                times = [0.09, 0.98, 1.21, 3.82, 2.77, 1.53, 1.57, 1.35, 1.51, 0.77]
                init_normal = [mean(times), std(times)]
                nll_trunc_normal(;params = init_weibull, times = times, trunc = [0, 4])

"""
function nll_trunc_normal(; params, times, trunc)
    mu, sigma = params
    if mu <= 0 || sigma <= 0
        return 1e10 #Inf
    end
    
    # Define distribution
    d0 = Normal(mu, sigma)
    #norm = cdf(d0, upper) # Only upper truncation
    #ll = sum(logpdf(d0, t) for t in times) - length(times) * log(norm)
    
    # Define truncated distribution
    lower, upper = trunc
    d = truncated( d0, lower, upper)
    
    # Remove zero times as they cause logpdf.(d, times) -> Inf || -Inf
    times_no_zeros = sort( filter(x -> x != 0, times) )
    
    # Compute log-likelihood
    ll = sum( sort(logpdf.(d, times_no_zeros)) ) #ll = sum( sort(logpdf.(d, times)) ) 
          
    return isfinite(-ll) ? -ll : 1e10
end


"""
Function:       dist_fit_plot

Description:    Plots data and fitted distribution

Arguments:      d0                  Fitted distribution, e.g. Gamma( shape = 2, scale = 1 )
                lower               Lower limit for truncation
                upper               Upper limit for truncation
                data_to_fit         Data that distribution is fitted to
                data_type::String   Description of data to be used in the plot title

Returns:    Plot with data as a histogram and line representing the fitted statistical distribution

Examples:   
            times = [0.09, 0.98, 1.21, 3.82, 2.77, 1.53, 1.57, 1.35, 1.51, 0.77]
            dist_fit_plot(; d0 = Gamma( 2, 1)
                          , lower = 0.0, upper = 60.0
                          , data_to_fit = times
                          , data_type = "tinf"
                          )
"""
function dist_fit_plot(; d0 
                        , lower , upper 
                        , data_to_fit
                        , data_type
                        )
    # Convert distribution to truncated version
    d = truncated( d0, lower, upper)
    # sample data points from fitted truncated distribution
    #d0_sample=rand(d0,10000)
    d_sample=rand(d,1000000)
        
    # Plot histogram of sampled data
    Plots.histogram( d_sample, alpha=0.5, normalize=:pdf, label="samples drawn from fitted dist")
    #histogram!(d0_sample,alpha=0.5, normalize=:pdf, color=:yellow)
    
    # Add simulated outbreak data to histogram
    Plots.histogram!(data_to_fit, normalize=:pdf, color=:red, alpha=0.5, label=data_type*" data")#" tinf or (treport-tinf) data")
    
    # Add simulated outbreak data with zero values removed to histogram
    data_to_fit_no_zeros = sort( filter(x -> x != 0, data_to_fit) )
    Plots.histogram!(data_to_fit_no_zeros, normalize=:pdf, color=:green, alpha=0.5, label=data_type*" data")#" tinf or (treport-tinf) data")
    
    # Overlay fitted (truncated) distribution PDF
    xs = range(0, upper; length=300)
    norm = cdf(d, upper) # normalization for truncation
    t_pdf_vals = pdf.(d, xs) ./ norm # truncated PDF
    Plots.plot!(xs, t_pdf_vals; label=data_type*" truncated fit", color=:blue, lw=2)
    
end


"""
Function: fit_multi_dist

Description:    Function to compare results for fitting Gamma and Weibull statistical distributions to data
                Note that this could potentially be extended to include other statistical distributions.

Arguments:  times           Vector of time values to be fit with statistical distribution
            init_gamma      Initial parameters (shape and scale as a two element vector) for fitting Gamma distribution
            init_weibull    Initial parameters (alpha and theta as a two element vector) for fitting Weibull distribution
            lower           Lower limit for truncation
            upper           Upper limit for truncation

Returns:    'fit_results' vector containing a dictionary. The Keys are the statistical distributions fitted to the data (times)
            and the values are the goodness of fit statistics (negative log-likelihood and AIC) and the optimised fit parameters
            for the statistical distributions.
            Example:
            1-element Vector{Dict{Any, Any}}:
            #    Distribution  =>  (Negative log-likelihood,    AIC,   [      Parameter 1, Parameter 2        ])
            Dict( "Gamma"   => (45.868899181683545, 95.73779836336709, [344.4013308355461, 0.26090968646639895])
                , "Weibull" => (1.0e10            , 2.0000000004e10  , [869.4059273903558, 0.09902739018403779])
                )

Examples:   
            #sims_simid_tinf_df = load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621/covidlike-1.3.14-sims_simid_tinf_nrep1200_1101621.jld2", "sims_simid_tinf_df")
            #times = sims_simid_tinf_df[1][:,2]
            #times2 = sample(times, 20) println(times2)
            times = [88.1, 85.9, 89.8, 85.5, 86.5, 83.7, 80.5, 83.8, 89.1, 85.4, 86.1, 89.2, 86.6, 89.1, 88.7, 86.5, 86.9, 78.8, 88.3, 83.4]
            init_gamma_1 = [mean(times)^2 / var(times), var(times)/mean(times)] #[2,1]
            init_weibull_1 = init_gamma_1

            fit_results = fit_multi_dist(; times = times2
                                         , init_gamma = init_gamma_1
                                         , init_weibull = init_weibull_1
                                         , lower = 0, upper = 90
                                         ) 
            # Returns
            1-element Vector{Dict{Any, Any}}:
            Dict( "Gamma"   => (45.868899181683545, 95.73779836336709, [344.4013308355461, 0.26090968646639895])
                , "Weibull" => (1.0e10            , 2.0000000004e10  , [869.4059273903558, 0.09902739018403779]))


"""
function fit_multi_dist(; times, init_gamma, init_weibull #,init_nbinom, init_lognorm, init_normal
                        , lower, upper) 
    
    # Fit all distributions and record results in a dictionary
    fit_results = Dict()

    # Gamma
    obj_f_gamma = init_gamma -> nll_trunc_gamma(;params = init_gamma, times = times, trunc = [lower,upper]) # [0, 60])
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


"""
Function:    erlang_truncated_means

Description     Splits an Erlang distribution into two parts and returns the mean x value of each part

Arguments   k::Int          Shape parameter for Erlang distribution (integer value makes it a special case of the Gamma distribution)  
            rate::Float64   Rate parameter for Erlang distribution
            lower_limit     Lower limit in x of the left hand split part of the distribution
            mid_limit       Upper limit in x of the left hand split part of the distribution and lower limit of the right hand split part
            upper_limit     Upper limit in x of the right hand split part of the distribution
Returns     Two Float64 values representing the mean x values of the two split parts of the distribution

Example     # Used to split the Erlang distribution reported in Knock et al (2021) fitted to the recovery time
            # for patients hospitalised with COVID-19. The mean of the whole distribution is 10.7 days. 
            # However, we want to split this into two: patients discharged within 24 hours and patients discharged
            # beyond 24 hours.
            erlang_truncated_means( k = 1, rate = 0.09, lower_limit = 0.0, mid_limit = 1.0, upper_limit = Inf)   
            # Returns (0.49250101230477183, 12.111111111110802)
            # Check also returns the mean from the original single distribution
            erlang_truncated_means( k = 1, rate = 0.09, lower_limit = 0.0, mid_limit = 10000, upper_limit = Inf)   
            # returns (11.11111, NaN)
            # Not quite the same as the 10.7 reported, but it is the same as E[t] = k / rate = 1/0.09 = 11.1111, so 
            # it is likely due to rounding of the rate paramter 1/0.0935 = 10.695

"""
function erlang_truncated_means(; k::Int
                                , rate::Float64
                                , lower_limit#::Union{Float64,Int,Inf}
                                , mid_limit#::Union{Float64,Int,Inf}
                                , upper_limit#::Union{Float64,Int,Inf}
                                )
    d = Gamma(k, 1/rate) # Erlang(k, λ)
    # pdf and cdf
    f(x) = pdf(d, x)
    F(x) = cdf(d, x)

    # (a) E[X | 0 ≤ X ≤ 1]
    num_a, _ = quadgk(x -> x * f(x), lower_limit, mid_limit) #quadgk(x -> x * f(x), 0.0, 1.0)
    den_a = F(mid_limit) - F(lower_limit)        # F(0.0) is 0 for Erlang
    mean_a = num_a / den_a

    # (b) E[X | X ≥ 1]
    num_b, _ = quadgk(x -> x * f(x), mid_limit, upper_limit) #quadgk(x -> x * f(x), 1.0, Inf)
    den_b = F(upper_limit) - F(mid_limit)
    mean_b = num_b / den_b

    return mean_a, mean_b
end


"""
Function    gamma_params_from_mode_cdf(mode_val, cdf_at_2; lower_shape=1.0+1e-6, upper_shape=1e3)

Description     Construct a `Gamma` distribution whose **mode** equals `mode_val` and whose
                cumulative distribution function (CDF) at `x = 2` matches `cdf_at_2` as closely
                as possible.

                This function solves for the Gamma shape parameter using a one-dimensional
                optimization over the domain `shape > 1`, which is required for the Gamma distribution
                to have a well-defined mode. For each candidate shape value, the corresponding scale

                mode = (shape - 1) * scale    ⇒    scale = mode_val / (shape - 1)

                The objective minimized is:

                abs( cdf(Gamma(shape, scale), 2) - cdf_at_2 )
                
Arguments   `mode_val::Real`:    Desired mode of the resulting Gamma distribution.
            `cdf_at_2::Real`:    Target value for the CDF at `x = 2`.
            `lower_shape::Real`: Lower bound for the shape parameter search (must be > 1).
            `upper_shape::Real`: Upper bound for the shape parameter search.

Returns   A `Gamma(shape, scale)` object from `Distributions.jl` whose mode equals
          `mode_val` and whose CDF at 2 is close to `cdf_at_2`.

Notes
        - The optimization may fail or converge to a boundary if the provided `mode_val`
          and `cdf_at_2` are incompatible with any Gamma distribution.
        - The returned distribution satisfies the mode constraint exactly (up to
          floating-point precision) and approximates the CDF constraint as closely as
          possible within the search bounds.

Example
        swab_time_gamma_d = gamma_params_from_mode_cdf(mode_val = 0.25, cdf_at_2 = 0.9
                                                        , lower_shape = 1.0 + 1e-6, upper_shape = 10)

        # Reproducibility
        Random.seed!(1234)

        # Draw samples
        N = 100_000
        x = rand(swab_time_gamma_d, N)

        # Quick diagnostics
        println("Sample mean ≈ ", mean(x))
        println("Sample variance ≈ ", var(x))
        println("Analytic mean = ", swab_time_gamma_d.α*swab_time_gamma_d.θ)
        println("Analytic variance = ", swab_time_gamma_d.α*swab_time_gamma_d.θ^2)
        println("Analytic mode = ", (swab_time_gamma_d.α - 1) * swab_time_gamma_d.θ)
        
        # Visualise
        Plots.histogram(x; bins=120, normalize=:pdf, alpha=0.5, label="ARI swab samples"
          , title="Gamma(swab_time_shape ≈ dollarsign(round(swab_time_gamma_d.α, digits=2)),\n swab_time_scale ≈ dollarsign(round(swab_time_gamma_d.θ,digits=2)))"
          , xlabel = "Time since hospital arrival (days)")
        xs = range(0, stop=6, length=600)
        plot!(xs, pdf.(swab_time_gamma_d, xs); label="Theoretical distribution", lw=2)
        vline!([0.25]; label="Mode ≈ 0.25", color=:red, lw=2, ls=:dash)


"""
function gamma_params_from_mode_cdf(;mode_val::Real, cdf_at_2::Real, lower_shape=1.0 + 1e-6, upper_shape=1e3)
    
    # --- Input validation ---
    @assert mode_val > 0 "mode_val must be > 0 so that scale = mode_val/(α-1) > 0."
    @assert 0.0 <= cdf_at_2 <= 1.0 "cdf_at_2 must be in [0, 1]."
    @assert lower_shape > 1.0 "lower_shape must be > 1."
    @assert upper_shape > lower_shape "upper_shape must exceed lower_shape."

    # Objective over shape > 1. Scale is determined by the mode constraint.
    obj(shape) = begin
        scale = mode_val / (shape - 1) # ensures given mode
        abs(cdf(Gamma(shape, scale), 2.0) - cdf_at_2)
    end

    res = optimize(obj, lower_shape, upper_shape, Brent())
    shape = Optim.minimizer(res)
    scale = mode_val / (shape - 1)
    return Gamma(shape, scale)
end
