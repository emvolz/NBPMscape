### Computing the additional care pathways for modelling HARISS sampling

## As per Knock et al (2021), Science Translational Medicine, 13(602)
# Hospitalised on general ward leading to recovery duration fitted with an 
# Erlang distribution Erlang(k=1, gamma=0.09), mean = k / gamma = 10.7 days (95% CI: 0.3, 39.4)
# We assume this represents COVID-19 patients admitted via ED attendance but excludes
# those that attended ED and were discharged instead of being admitted.
# Therefore, the proportion of the fitted Erlang distribution with a duration <24h (<1 day)
# represents the admissions categorised as 'short stay'.
# Saigal et al. (2025), BMC Emergency Medicine, 25:11, reported a cohort of COVID-19 patients 
# attending ED with 66% discharged within 24h and 34% with a longer stay.
# The 66% includes those discharged without admission. If the proportion of the Erlang
# distribution < 1 day = X%, then the 34% is proportional to the 1-X% of those admitted 
# to hospital in Knock et al (2021). We can then compute a normalisation factor from this
# and apply this to the 66% that were discharged. From the current model, using analysis 
# of simulation results, we know that ~3.4% of infections are admitted to hospital but
# not ICU.

## Compute proportion of hospital admissions with recovery periods (assume this is equivalent 
# to discharge) < 24h (1-day)

using Distributions
using Plots

# Parameters (Erlang distribution is a special case of the Gamma distribution
# with integer shape parameter value)
k     = 1          # shape from Knock et al (2021)
rate  = 0.09       # rate from Knock et al (2021)
x0    = 1.0        # Split distribution at x = 1day (or 24h)

dist = Gamma(k, 1/rate)   # Gamma(shape=k, scale=1/rate) ⇒ Erlang(k, rate)

# Proportion of the distribution below x = 1day (24h)
p_admitted_but_discharged_within_24h = cdf(dist, x0) # 0.086 = 8.6% of all COVID-19 hospital admissions (assume does not include those discharged directly from Emergency Department ED)
println("P(X ≤ $x0) = $p_below")

# Split between short stay and long stay for COVID-19 patients from Saigal et al (2025)
p_ED_discharged_within_24h = 0.66
p_ED_discharged_after_24h = 0.34

# Need a normalisation factor because trying to subtract the 'short stay' 
# admissions from the total discharged within 24h, the latter includes
# those who were discharged and not admitted
norm_factor = (1-p_admitted_but_discharged_within_24h)/ p_ED_discharged_after_24h

# Proportion of COVID-19 positive that attend ED and are discharged directly
# (i.e. without admission to another ward)
# prop discharged from ED = prop discharged within 24h - prop admitted AND discharged within 24h
p_ED_discharged_directly = p_ED_discharged_within_24h - (norm_factor*p_admitted_but_discharged_within_24h)

# Of those infected with COVID-19 that visit hospital (assume via ED first):
# ~43% discharged directly from ED
# ~23% admitted to ward but discharged within 24h
# ~34% admitted to ward and discharged later than 24h
prop_severe_hosp_short_stay = 23 / (23+34)
prop_severe_hosp_short_stay = 34 / (23+34)

# ~3.4% of all COVID-19-like infections in model simulation results are admitted to hospital
# Considering the breakdown above, this is equivalent to 
# 1 - p_ED_discharged_directly = 0.571 = 57.1% of those attending ED
p_attend_ED_but_discharged = (0.034 / (1 - p_ED_discharged_directly)) - 0.034
# Therefore, ~2.6% of all COVID-19-like infected individuals attend ED but are discharged
# We take this from the weighting previously allocated to mild symptoms,
# which was a balancing weight.

# Distribution of severity from previous simulation (median values from 1000 simulations),
# as computed using inf_severity_estimate function in
# 'src/estimate_severity_age_weights.jl', was:
# - Asymptomatic 52.3%
# - Mild 39.0%
# - Moderate (GP visit) 4.7%
# - Severe (hospital admission but not ICU) 3.4%
# - Very severe (ICU - possibly via GP and hospital) 0.5%
# So as the allocation of infection severity is set up, which allocates
# ICU and hospital admission based on age disaggregated probabilities and
# allocates the remainder using non-age disaggregated probability 
# estimates, means that the 2.6% of total cases should be converted
# to the 0.026 / (1 - asymptomatic - severe_wt - very_severe_wt)
# Based on weights computing from previous simulation using model with just 5 categories (asymptomatic, mild, moderate, severe, very severe)
# Symptomatic = 1 - 0.523 = 0.477
# So,
## Proportion of all infections
# moderate_ED (proportion of all infections) = 0.026
# moderate_GP (proportion of all infections) = (1-0.523) * 0.10825 = 0.05164 ~ 4.7% (as above) # 0.111 is computed from FluSurvey of people with ILI symptoms. GP visits (in-person and phone) / non-hospital options = (9+1.825)/(100-2.8)
# mild (proportion of all infections) = 1 - 0.523 - 0.052947 - 0.026 - 0.034 - 0.005 = 0.35905 ~ 36.4% = 39.0% - 2.6% (as above)
## Proportion of symptomatic infections
# moderate_ED (proportion of symptomatic) = 0.026 / ( 1 - 0.523 ) = 0.05451
# moderate_GP (proportion of symptomatic) = 0.10825 # computed from FluSurvey. GP visits (in-person and phone) 10.825% = (9% + 1.825%)
# mild (proportion of symptomatic) = 0.35905 / ( 1 - 0.523 ) = 0.75273
## Proportion of symptomatic infections after removing severe and very severe infections
# moderate_ED (proportion of symptomatic after removing asymptomatc, severe and very severe) = 0.05451 / ( 1 - 0.523 - 0.034 - 0.005) = 0.12445
# moderate_GP (proportion of symptomatic after removing asymptomatc, severe and very severe) = 0.05164 / ( 1 - 0.523 - 0.034 - 0.005) = 0.11790
# mild (proportion of symptomatic after removing asymptomatc, severe and very severe) = 1 - 0.12445 - 0.11790 = 0.75765 # This is the remainder


## Plot Erlang distribution split in two - short stay and long stay hospital admissions

# Prepare grid for plotting
xmax = 80.0                  # upper limit for x-axis
xs   = range(0, xmax; length=500)
ys   = pdf.(dist, xs)

# Split x-range into below/above x0
xs_below = xs[xs .<= x0]
ys_below = pdf.(dist, xs_below)

xs_above = xs[xs .>= x0]
ys_above = pdf.(dist, xs_above)

# Base plot of the PDF
plt = Plots.plot(xs, ys,
           label = "PDF",
           xlabel = "x",
           ylabel = "f(x)",
           linewidth = 2,
           legend = :topright)

# Shade area below x0
Plots.plot!(plt, xs_below, ys_below,
      fillrange = 0,
      fillalpha = 0.4,
      color = :blue,
      label = "Area ≤ $x0")

# Shade area above x0
Plots.plot!(plt, xs_above, ys_above,
      fillrange = 0,
      fillalpha = 0.4,
      color = :red,
      label = "Area > $x0")

# Vertical line at x0
vline!(plt, [x0],
       color = :black,
       linestyle = :dash,
       label = "x = $x0")

display(plt)