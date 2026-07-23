#=  
Functions for comparing results from different sampling strategies for 
time to detection (TD) and cumulative infections at TD

There are two methods, both of which seek to select the strategy with the minimum value
for a function of the format: mean + penalty x risk measure
The two methods used are:
- Conditional value-at-risk (CVaR) function: E[x] + phi    * CVaRα[x]
- Loss function:                             E[x] + lambda * variance[x]

This file includes the following functions:
- mean_risk_loss    Computes the mean-risk objective, J = mean(x) + (lambda * variance(x))
                    and returns a named tuple including this as well as some input and intermediate values
- mean_risk_cvar    Computes the mean-risk objective using Conditional Value-at-Risk (CVaR).
                    Uses the Rockafellar-Uryasev (2002) formulation as also used in Yin et al (2023):
                    CVaRα(z) = min_{η} { η + 1/(1-α) * E([z - η]⁺) } 
                             = min_{η} { η + 1/(1-α) * E(max([z - η],0) }
                    The mean-risk objective is: C = E(z) + ϕ * CVaRα(z)
                    Returns a named tuple including the mean-risk objective as well as some input and intermediate values
- strategy_rank     Uses the {mean_risk_loss} and {mean_risk_cvar} functions to rank
                    multiple sampling strategies based on the selected mean-risk objective
                    and user input parameter values

And the following helper functions for preparing inputs:
- 

=#

"""
# Function      mean_risk_loss(; x::Vector{<:Real}, lambda::Real)

# Description   Computes the mean-risk objective, J:
                J = mean(x) + (lambda * variance(x))
                Note that the variance risk measure is two-tailed and therefore requires the
                distribution of x to be approximately normal and symmetric

# Arguments     - x::Vector{<:Real}` : Vector of simulated values (i.e. times to detection (TD))
                - lambdaϕ::Real      : Risk weight coefficient to penalise risk
                                        
# Returns       A `NamedTuple` with fields:
                - `expected_value (mean)`   : mean(x)
                - `variance`                : variance(x)
                - 'risk penalty'            : lambda
                - `mean_risk_obj`   : Mean-risk objective = mean(x) + (variance(x) * lambda)

# Example
                x = [2.0, 4.0, 6.0, 9.0, 14.0]
                result_0 = mean_risk_loss(x = x, lambda = 0)
                # (mean_risk_obj = 7.0, mean_x = 7.0, variance_x = 22.0, lambda = 0)
                result_05 = mean_risk_loss(x = x, lambda = 0.5)
                # (mean_risk_obj = 18.0, mean_x = 7.0, variance_x = 22.0, lambda = 0.5)

"""
function mean_risk_loss(; x::Vector{<:Real}, lambda::Real)

    # ── Input Validation ────────────────────────────────────────────────────
    isempty(x)          && throw(ArgumentError("Input vector `x` must not be empty."))
    !(lambda >= 0)       && throw(ArgumentError("Risk weight `lambda` must be equal to or greater than 0. Got lambda = $lambda"))

    # ── Compute mean-risk objective ─────────────────────────────────────────
    mean_x = mean(x)
    variance_x = var(x)
    mean_risk_obj = mean_x + (variance_x * lambda)

    return (
          mean_risk_obj = mean_risk_obj
        , mean_x        = mean_x
        , variance_x    = variance_x
        , lambda        = lambda
    )
end


"""
# Function      mean_risk_cvar(; z::Vector{<:Real}, phi::Real, alpha::Real, VaR_quantile_method::String = "exact")

# Description   Compute the mean-risk objective using Conditional Value-at-Risk (CVaR).

                Uses the Rockafellar-Uryasev (2002) formulation as also used in Yin et al (2023):
                CVaRα(z) = min_{η} { η + 1/(1-α) * E([z - η]⁺) } 
                         = min_{η} { η + 1/(1-α) * E(max([z - η],0) }

                The mean-risk objective is:
                J = E(z) + ϕ * CVaRα(z)

# Arguments     - `z::Vector{<:Real}`       : Vector of observed/simulated durations or costs
                - `ϕ::Real`                 : Risk weight coefficient ∈ [0, 1], where:
                                                0 = risk-neutral (expected value only)
                                                1 = fully risk-averse (maximum tail penalty)
                - `α::Real`                 : Confidence level, must be ∈ [0, 1]
                - 'VaR_quantile_method::String' : ["exact","interpolate"] Choice between returning an exact 
                                                  value within the vector z or potentially interpolating between 
                                                  two values. The "exact" option is closer to the definition of
                                                  VaR used in Yin et al (2023) and Rockafellar & Uryasev (2002).
                                                  Also see Hyndman & Fan (1996) for further details on quantile
                                                  calculation.

# References    Hyndman & Fan (1996), Sample Quantiles in Statistical Packages, The American Statistician, 50(4):361-365
                Rockafellar & Uryasev (2002), Conditional value-at-risk for general loss distributions, Journal of Banking & Finance, 26(7):1443-1471 
                Yin et  al (2023), COVID-19: Data-Driven optimal allocation of ventilator supply under uncertainty and risk, European Journal of Operational Research, 304(1):255-275

# Returns       A `NamedTuple` with fields:
                - `mean_risk_obj`       : Mean-risk objective = E(z) + ϕ * CVaRα(z)
                - `expected_value`      : E(z) = mean(z)
                - 'penalty_factor'      : phi
                - 'alpha_quantile'      : alpha
                - `VaR`                 : Value-at-Risk at level α
                - `CVaR`                : Conditional Value-at-Risk at level α
                - 'VaR_quantile_method' : VaR_quantile_method
                
# Example
                x = [2.0, 4.0, 6.0, 9.0, 14.0]
                result_05_80 = mean_risk_cvar(z = x, phi = 0.5, alpha = 0.80)
                # (mean_risk_obj = 14.0, expected_value = 7.0, VaR = 10.0, CVaR = 14.0, phi = 0.5, alpha = 0.8)
                result_05_95 = mean_risk_cvar(z = x, phi = 0.5, alpha = 0.95)
                # (mean_risk_obj = 15.499999999999998, expected_value = 7.0, VaR = 13.0, CVaR = 16.999999999999996, phi = 0.5, alpha = 0.95)
"""
function mean_risk_cvar(; z::Vector{<:Real}, phi::Real, alpha::Real, var_quantile_method::String = "exact") #VaR_quantile_method = "interpolate"

    # ── Input Validation ────────────────────────────────────────────────────
    isempty(z)          && throw(ArgumentError("Input vector `z` must not be empty."))
    !(0 < alpha < 1)    && throw(ArgumentError("Confidence level `alpha` must be in [0, 1]. Got alpha = $alpha"))
    !(0 ≤ phi ≤ 1)      && throw(ArgumentError("Risk weight `phi` must be in [0, 1]. Got phi = $phi"))

    n = length(z)

    # ── Step 1: Expected Value ───────────────────────────────────────────────
    expected_val = mean(z)

    # ── Step 2: VaR — the α-quantile of z ───────────────────────────────────
    # quantile() uses linear interpolation by default in Julia (type=7)
    # Rockafellar-Uryasev formulation uses the empirical CDF so we set beta=1:
    #   VaRα = inf{η ∈ R : Fz(η) ≥ α}
    VaR = quantile(  z  # Vector of values to return VaR (quantile) for
                   , alpha # quantile to return
                   , sorted=false # Does not assume that vector is already sorted (this is also the default)
                   , alpha = var_quantile_method == "exact" ? 0 : 1 # this value in the quantile function determines whether an exact or interpolated value are returned from the quantile function
                   , beta = 1  # second parameter for controlling quantile calculation
                   )

    # ── Step 3: CVaR via Rockafellar-Uryasev Formula ────────────────────────
    #   CVaRα(z) = η + 1/(1-α) * E([z - η]⁺)
    #   where η = VaR and [a]⁺ = max(a, 0)
    exceedances   = max.(z .- VaR, 0)           # [z - VaR]⁺ for each scenario
    expected_exc  = mean(exceedances)            # E([z - VaR]⁺)
    CVaR          = VaR + (1 / (1 - alpha)) * expected_exc

    # ── Step 4: Mean-Risk Objective ──────────────────────────────────────────
    #   J = E(z) + ϕ * CVaRα(z)
    mean_risk_obj = expected_val + phi * CVaR

    return (
        mean_risk_obj  = mean_risk_obj    
        ,expected_value = expected_val
        ,VaR            = VaR
        ,CVaR           = CVaR
        ,phi            = phi
        ,alpha          = alpha
    )
end


"""
# Function      strategy_rank(; df::DataFrame, rank_by::String, lambda::Float64, phi::Float64, alpha::Float64)

# Description   Compute mean-risk objectives using two methods:
                (1) Conditional Value-at-Risk (CVaR) via {mean_risk_cvar} function
                (2) Loss function using variance as risk measure via {mean_risk_loss} function
                for one or more vectors of data organised in a dataframe.
                The output is a df containing labels for the strategies (taken from the input df column names)
                in the first column and mean-risk objective data in the other columns.
                The strategies will also be ordered (ranked) based on the strategy selected 
                in the 'rank_by' argument.
                A plot will also be produced showing the "rank_by" risk measure against the
                expected values for each strategy.
                Note: that the variance risk measure used in the (2) Loss function is two-tailed and
                therefore requires the distribution of x to be approximately normal and symmetric
                    
# Arguments     - df::DataFrame     : DataFrame containing the data to be analysed in columns with the strategy labels as the column names
                - rank_by::String   : Options are ["cvar","loss"] which are the two different mean-risk objective values that can be calculated
                - lambda::Float64   : risk penalty factor in the mean_risk_loss function
                - phi::Float64      : risk penalty factor in the mean_risk_cvar function
                - alpha::Float64    : Quantile used in the conditional value-at-risk calculation (via mean_risk_cvar function)
                - var_quantile_method::String : ["exact","interpolate"] Choice between returning an exact 
                                                  value within the vector z or potentially interpolating between 
                                                  two values. The "exact" option is closer to the definition of
                                                  VaR used in Yin et al (2023) and Rockafellar & Uryasev (2002).
                                                  Also see Hyndman & Fan (1996) for further details on quantile
                                                  calculation.

# Returns       A dataframe, sorted by the 'rank_by' mean-risk objective, with columns:
                - `strategy_name`                       : Name of strategy using the column names of the input df
                - `expected_value`                      : E(z) = mean(z)
                - 'variance'                            : variance of vector
                - 'penalty_factor_loss_function'        : phi
                - 'mean_risk_objective_loss_function'   : J(theta) = E[x] + penalty_factor_loss_function * variance(x)
                - 'penalty_factor_cvar'                 : phi
                - 'alpha_quantile_cvar'                 : alpha
                - `VaR`                                 : Value-at-Risk at level α
                - `CVaR`                                : Conditional Value-at-Risk at level α
                - 'mean_risk_objective_cvar'            : C(θ) = E[x|θ] + phi * CVaRα(x|θ)
                
# Examples
                input_df_1 = DataFrame( strategy_a = [2.0, 4.0, 6.0, 9.0, 14.0]
                                      , strategy_b = [12.0, 14.0, 16.0, 19.0, 24.0]
                                      , strategy_c = [0.2, 0.4, 0.6, 0.9, 1.4]
                                      )
                result_1 = strategy_rank( df = input_df_1, rank_by = "cvar"
                                        , lambda = 0.5, phi = 0.5
                                        , alpha = 0.80, var_quantile_method = "exact")
                
                input_df_2 = DataFrame( strategy_norm_a = max.( rand( Normal(40, 10), 1000), 0 )
                                      , strategy_norm_b = max.( rand( Normal(40, 20), 1000), 0 )
                                      , strategy_norm_c = max.( rand( Normal(50, 5), 1000), 0 )
                                      )
                result_2 = strategy_rank( df = input_df_2, rank_by = "cvar"
                                        , lambda = 0.5, phi = 0.5
                                        , alpha = 0.80, var_quantile_method = "exact")
                
"""
function strategy_rank(; df::DataFrame, rank_by::String, lambda::Float64, phi::Float64, alpha::Float64, var_quantile_method::String = "exact")

    # ── Input Validation ────────────────────────────────────────────────────
    isempty(df)         && throw(ArgumentError("Input vector `z` must not be empty."))
    !(0 < alpha < 1)    && throw(ArgumentError("Confidence level `alpha` must be in [0, 1]. Got alpha = $alpha"))
    !(0 ≤ phi ≤ 1)      && throw(ArgumentError("Risk weight `phi` must be in [0, 1]. Got phi = $phi"))
    !(lambda >= 0)      && throw(ArgumentError("Risk weight `lambda` must be equal to or greater than 0. Got lambda = $lambda"))
    !(var_quantile_method in ["exact","interpolate"]) && throw(ArgumentError("VaR quantile method `var_quantile_method` must be either 'exact' or 'interpolate'. Got var_quantile_method = $var_quantile_method"))

    # Initialise output dataframe
    empty_col_vec = Vector{Union{Missing, Float64}}(missing, ncol(df))
    results_df = DataFrame( 
                              strategy_name = names(df)
                            , expected_value = copy(empty_col_vec)
                            , variance       = copy(empty_col_vec)
                            , penalty_factor_loss_function = copy(empty_col_vec)
                            , mean_risk_objective_loss_function = copy(empty_col_vec)
                            , penalty_factor_cvar = copy(empty_col_vec)
                            , alpha_quantile_cvar = copy(empty_col_vec)
                            , VaR = copy(empty_col_vec)
                            , CVaR = copy(empty_col_vec)
                            , mean_risk_objective_cvar = copy(empty_col_vec)
    )

    # Print warning for two-tailed risk measure
    println("**Note that the mean-risk objective for the loss function uses variance as the risk measure which is two-tailed and therefore requires the distribution to be approximately normal and symmetric**")
    
    # Loop through input df and compute mean-risk objective using both methods for each strategy (column)
    # and add results to output df
    for x in 1:size(df,2) #x=1
    
        # Compute mean-risk objectives
        strat_data = df[:,x]
        mro_cvar = mean_risk_cvar( z = strat_data, phi = phi, alpha = alpha, var_quantile_method = var_quantile_method)
        mro_loss = mean_risk_loss( x = strat_data, lambda = lambda)

        # Add results to output df
        results_df[x,"expected_value"]                    = mro_cvar.expected_value # This should be exactly the same as in mro_cvar.expected_value and mro_loss.mean_x
        results_df[x,"variance"]                          = mro_loss.variance_x
        results_df[x,"penalty_factor_loss_function"]      = mro_loss.lambda
        results_df[x,"mean_risk_objective_loss_function"] = mro_loss.mean_risk_obj
        results_df[x,"penalty_factor_cvar"]               = mro_cvar.phi
        results_df[x,"alpha_quantile_cvar"]               = mro_cvar.alpha
        results_df[x,"VaR"]                               = mro_cvar.VaR
        results_df[x,"CVaR"]                              = mro_cvar.CVaR
        results_df[x,"mean_risk_objective_cvar"]          = mro_cvar.mean_risk_obj

    end
        
    # Sort output df by the selected mean-risk objective for each strategy
    if     rank_by == "cvar"
        results_df_sorted = sort( results_df, :mean_risk_objective_cvar)
    elseif rank_by == "loss"
        results_df_sorted = sort( results_df, :mean_risk_objective_loss_function)
    end

    # Plot raw risk measures against expected values for each strategy
    p_raw = Plots.scatter( results_df.expected_value, results_df.variance, label="Variance used in loss function"
                    , xlabel = "Expected value", ylabel = "Raw risk measure")
    p_raw = Plots.scatter!(results_df.expected_value, results_df.CVaR, label="CVaR uesd in CVaR function", color = :black)

    # Mark in red the top ranked strategy in terms of mean-risk objective
    if rank_by == "cvar"
        p_adj = Plots.scatter!([results_df_sorted.expected_value[1]], [results_df_sorted.CVaR[1]]
                                , label="Top ranked strategy using the CVaR function", color = :red)
    elseif rank_by == "loss"
        p_adj = Plots.scatter!([results_df_sorted.expected_value[1]], [results_df_sorted.variance[1]]
                                , label="Top ranked strategy using loss function", color = :red)
    end

    
    # Plot penalised risk measures against expected values for each strategy
    p_adj = Plots.scatter( results_df.expected_value, results_df.mean_risk_objective_loss_function, label="mean-risk objective for loss function"
                    , xlabel = "Expected value (E[x])", ylabel = "mean-risk objective (=E[x] + penalty factor * risk measure)")
    p_adj = Plots.scatter!(results_df.expected_value, results_df.mean_risk_objective_cvar, label="mean-risk objective for CVaR function"
                            , color = :black)
    # Mark in red the top ranked strategy in terms of mean-risk objective
    if rank_by == "cvar"
        p_adj = Plots.scatter!([results_df_sorted.expected_value[1]], [results_df_sorted.mean_risk_objective_cvar[1]], label="Top ranked strategy CVaR function"
                                , color = :red)
    elseif rank_by == "loss"
        p_adj = Plots.scatter!([results_df_sorted.expected_value[1]], [results_df_sorted.mean_risk_objective_loss_function[1]], label="Top ranked strategy using loss function"
                                , color = :red)
    end

    p_combined = Plots.plot( p_raw, p_adj, layout = (1, 2)
                            , size = (800,500)) #default size is 600,400
    display( p_combined )
    
    println("Strategy rankings")
    println("_________________")
    println("$(show(results_df_sorted[:,["strategy_name","expected_value","mean_risk_objective_loss_function","mean_risk_objective_cvar"]]; allcols=true, truncate=0))")
        
    return ( results = (results_df_sorted = results_df_sorted, mean_risk_plot = p_combined)  )
end



# Helper function to convert dictionary data (formatted with strategy names as keys and
# dataframes containing sampling results as columns) to a df (formatted with strategy names as column names and
# the selected data type as row values) for input into strategy_rank function
using DataFrames
function results_dict_to_df(; df_dict::Dict, col_to_extract = :ICU_TD )
    # Initialise an empty results DataFrame
    results_df = DataFrame()

    # Loop through each key-value pair in the dictionary
    for (key, df) in df_dict
        # Extract the desired column and insert it into the results DataFrame
        # using the dictionary key as the new column name
        results_df[!, key] = df[!, col_to_extract]
    end
    
    return(results_df)
end

# Helper function to filter strategies (using strings contained in the name) before
# running the strategy_rank function
# Filter df for columns containing all substrings in the column name
function filter_strat_names(; df::DataFrame, substrings::Vector{String} )
    cols = filter(
        #name -> any(sub -> occursin(sub, String(name)), substrings),
        name -> all(sub -> occursin(sub, String(name)), substrings),
        names(df)
    )

    return( select(df, cols) )
end