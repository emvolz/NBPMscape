# NBPMscape.jl — Parameters
Updated: 30 July 2026
> **Disease Transmission and Cluster Growth Simulation**
>
> This document describes the parameters used in simulating
> disease transmission and cluster growth for a respiratory pathogen within England.
>
> The pathogen specific parameter values are for a SARS-CoV-2-like respiratory
> pathogen and are the default values used in NBPMscape.jl. The default parameters are defined in `src/core.jl` by `create_default_parameters()`.
>
> The pathogen specific and sampling default parameters can be overridden via a YAML configuration file (e.g.
> `config/outbreak_params_covid19_like.yaml`). However, parameter values controlling transmission dynamics are not currently updated from configuration files.
>
> **How to use:** Pathogen specific and sampling parameters can be set via a YAML config file and loaded with:
```julia
 initialize_parameters("config/outbreak_params_covid19_like.yaml") 
 ```
> Default parameters are returned by:
```julia
initialize_parameters()
```
> with no argument.

---

## 1. Importation

We assume that the initial case of infection is imported from outside the UK and
that imports continue after the initial case. 
There are three models available to select for generating imported infections. 
The first two follow a t-distribution with the parameters listed below. The difference is the level of truncation of the distribution. 'importmodel = :TDist_tlb_1' is truncated at the quantile '1/nimports=1/1000=0.001' while 'importmodel = :TDist_tlb_2' is truncated at the quantile = 0.01. The third option, 'importmodel = :GLEAM_outbreak_country_not_specified' takes Poisson samples from a daily time series of mean daily importations generated using the GLEAM model [St-Onge et al (2025) and Reddy et al (2026)].

The default parameters for the importation methods using a t-distribution assume a daily importation rate of 0.5 and 1,000 total imports, as per the analysis by du Plessis et al. (2020) of the establishment of the SARS-CoV-2 epidemic in the UK. t-distribution parameters are derived from a personal analysis of lineages
studied in Volz et al. (2020).

```julia
importrate  = 0.5   # Daily importation rate (constant rate model)
nimports    = 1000  # Total number of imports (du Plessis et al. 2020)
import_t_df = 2.48  # t-distribution degrees of freedom (Volz et al. 2020 lineage analysis)
import_t_s  = 8.75  # t-distribution scale parameter (Volz et al. 2020 lineage analysis)
initial_dow = 1     # Day of week of initial import (Sunday = 1, ..., Saturday = 7)
```
The location of the imported cases in the UK is determined using a probability based on international air passenger arrivals and airport location. International air passenger traffic is sourced from the [UK Civil Aviation Authority Table 12.1](https://www.caa.co.uk/data-and-analysis/uk-aviation-market/airports/uk-airport-data/) for 2024. [Accessed: 31 October 2024 and 15 March 2025]. Geographic area definitions from the UK Office for National Statistics (ONS) were used to determine the International Territory Level 2 (ITL2) in which airports in England are located.

---

## 2. Infectivity

Infectivity varies with time since becoming infectious and is modelled using a
Gamma distribution. The shape and scale parameters are manually calibrated to give a generation time of
5–6 days, consistent with Hart et al. (2022), Lau et al. (2021), Chen et al.
(2022), Xu et al. (2023), and Bi et al. (2020). The overall scale of
transmission is set by `infectivity`. This is manually tuned to give R ~ 2 using the 'infectivitytoR()' function to compute R for 
test simulations.

```julia
infectivity       = 2.00
infectivity_shape = 1.65
infectivity_scale = 1.875
```

---

## 3. Latent Period Duration

The latent infection period is the duration between the time of infection and the time that an individual becomes infectious. This differs from the incubation period which is the duration between the time of infection and the onset of symptoms. These two periods may be different, which can result in pre-symptomatic transmission. 

We assume that the latent infection period for an individual is drawn from a Gamma distribution. We use the Gamma distribution parameters below which we infer from the mean latent period estimate of 3.3 days (95% CI: 0.2–7.9) by Zhao et al (2020).

```julia
latent_shape = 3.26
latent_scale = 0.979
```

---

## 4. Infectious Period Duration

Individuals remain infectious for different periods of time. The length of the infectious period follows a Gamma distribution. We use the parameters estimated by Verity et al (2020) which have a mean of 24 days: Gamma(shape = 8.16, scale = 3.03).

```julia
infectious_shape = 8.16
infectious_scale = 3.03
```

---

## 5. Number of Contacts (Age-Stratified)

The number of contacts varies according to the individual's age and the contact setting. 
We analysed the POLYMOD UK survey [Massong et al, (2008)] to produce an empirical distribution of the number of contacts in various settings. These were fitted to statistical distributions using the *fitdistrplus* R package [Delignette-Muller & Dutang (2015)].
Contact matrices are generated using the POLYMOD UK survey (Mossong et al. 2008) via the R package *socialmixr*.
We use three contact settings: home, work/school and other, which comprises the travel, leisure and other categories from the original survey.

All distribution parameters are stored as **age-indexed vectors** (element 1 =
age 0 years, up to age 100), loaded from `CONTACT_DISTRIBUTIONS` data at runtime.

### 5.1 Household contacts (F-links)

Negative Binomial distribution, parameters vary by age (POLYMOD home setting):

```julia
fnegbinomr # Vector: NegBin r parameter for household contacts, by single year age
fnegbinomp # Vector: NegBin p parameter for household contacts, by single year age
```

### 5.2 Workplace / school contacts (G-links)

Negative Binomial distribution, parameters vary by age (POLYMOD work/school setting):

```julia
gnegbinomr # Vector: NegBin r parameter for work/school contacts, by single year age
gnegbinomp # Vector: NegBin p parameter for work/school contacts, by single year age
```

### 5.3 Casual / other contacts (H-links)

Gamma distribution, parameters vary by age (POLYMOD combined transport, leisure,
other settings):

```julia
oorateshape  # Vector: Gamma shape for other contacts, by single year age
ooratescale  # Vector: Gamma scale for other contacts, by single year age
oorateshape1 # Vector: = oorateshape + 1 (excess rate distribution)
ooratescale1 # Vector: = ooratescale (same scale, used with oorateshape1)
```

### 5.4 Age-assortative contact matrices

Single-year-age contact matrices are derived from POLYMOD UK [Mossong et al. (2008)] for each setting and used to sample the age of secondary cases conditional on the age of the infector:

```julia
f_contact_matrix_age # Home contact matrix (101 × 101, single year ages 0–100)
g_contact_matrix_age # Work/school contact matrix (101 × 101)
o_contact_matrix_age # Other contact matrix (101 × 101)
```

---

## 6. Changes to Contact Networks

The contact network has a dynamic structure and allows for contacts to be made and lost. 
The rate at which contacts are made and lost varies by the type of contact.:

```julia
frate = 0.0      # Rate of gaining/losing household (F) contacts
grate = 1 / 30.0 # Rate of gaining/losing work/school (G) contacts
```

---

## 7. Commuting

The model incorporates movement of individuals between regions. The commuting rates are based on an analysis of 'origin-destination' data for England and Wales which is part of the ONS UK 2021 census and available through the ([nomis website](https://www.nomisweb.co.uk/sources/census_2021_od) [Accessed: 14 March 2025]). In particular, the 'ODWP01EW - Location of usual residence and place of work' dataset was used. Commuting applies only to individuals aged 16–64.

```julia
commuterate = 2.0 # Commuting rate parameter
```

---

## 8. Relative Probability of Transmission by Contact Type

Relative transmission probabilities by contact type are normalised to household
contact hours (10.76 hours) from Danon et al. (2013), Table S2:

```julia
fcont  = 10.76 / 10.76 # ≈ 1.000 — household contacts (F-links)
gcont  =  6.71 / 10.76 # ≈ 0.623 — work/school contacts (G-links)
oocont =  8.09 / 10.76 # ≈ 0.751 — casual/other contacts (H-links)
```

---

## 9. Day-of-Week Variation in Transmission

Contact rates are scaled by the day of the week, sourced from the POLYMOD UK
survey (Mossong et al. 2008):

```julia
#        Sun        Mon        Tue        Wed        Thu        Fri        Sat
dowcont = (0.1043502, 0.1402675, 0.1735913, 0.1437642, 0.1596205, 0.1445298, 0.1338766)
```

---

## 10. Transmission Reduction

### 10.1 Hospitalised individuals

Transmission is reduced for individuals admitted to hospital (general ward), ICU,
or stepdown ward, due to reduced contacts and use of protective equipment:

```julia
ρ_hosp = 0.250 # Transmission rate is 25% of normal while hospitalised/in ICU/stepdown
```

### 10.2 Asymptomatic individuals

Asymptomatic individuals have a reduced transmission rate, consistent with value for COVID-19 in Knock
et al. (2021) Supplementary Information:

```julia
ρ_asymptomatic = 0.223 # Transmission rate is 22.3% of normal for asymptomatic individuals
```

---

## 11. Infection Severity (Age-Stratified)

Infection severity is determined probabilistically, stratified by the age of the
infectee. Probabilities are sourced from Knock et al. (2021) and Saigal et al.
(2025). The possible severity categories are:

| Severity                    | Description                                                  |
|-----------------------------|--------------------------------------------------------------|
| `:asymptomatic`             | No symptoms; reduced transmission rate                       |
| `:mild`                     | Symptomatic; no healthcare contact                           |
| `:moderate_GP`              | Symptomatic; visits GP only                                  |
| `:moderate_ED`              | Symptomatic; visits Emergency Department but not admitted    |
| `:severe_hosp_short_stay`   | Admitted to hospital; stay < 24 hours                        |
| `:severe_hosp_long_stay`    | Admitted to hospital; stay ≥ 24 hours                        |
| `:verysevere`               | Admitted to ICU                                              |

Age-stratified probabilities are read from data files at runtime:

```julia
symptomatic_IHR_IFR_by_age_file = "data/severity_care_probabilities/covid_knock_symptomatic_ihr_ifr_by_age.csv"
care_pathway_prob_by_age_file   = "data/severity_care_probabilities/covid_knock_care_pathway_prob_by_age.csv"
```
These files can also be changed in the configuration files described at the beginning.


Among symptomatic non-ICU hospital admissions, the split between short and long
stay is estimated based on Knock et al. (2021) and Saigal et al. (2025):

```julia
prop_severe_hosp_short_stay = 0.40351 # < 24h stay
prop_severe_hosp_long_stay  = 0.59649 # ≥ 24h stay
```

Among symptomatic individuals not admitted to hospital:

```julia
prop_moderate_ED = 0.12445 # Visits ED but discharged without admission. Estimated based on Knock et al. (2021) and Saigal et al. (2025). 
prop_moderate_GP = 0.11790 # Visits GP only (UKHSA (2025b) FluSurvey 2024/25)
prop_mild        = 0.75765 # = 1 - prop_moderate_GP - prop_moderate_ED
```

---

## 12. Care Pathway Rates

Rates (= 1 / mean duration) for progression through the care pathway. Sources:
Knock et al. (2021), Docherty et al. (2020), Saigal et al. (2025), and UKHSA (2025a) Data quality report: national flu and COVID-19 surveillance report (27 May 2025).

### 12.1 Rates from symptom onset to first care contact

```julia
gp_only_rate            = 1 / 5 # Lag from ARI symptom onset to GP visit (3–7 days typical). UKHSA (2025a) Data quality report: national flu and COVID-19 surveillance report (27 May 2025).
ed_direct_rate          = 1 / 5 # Lag from symptom onset to ED visit (direct, no GP)
gp_before_hosp_rate     = 1 / 3 # Lag from symptom onset to GP visit (before onward referral)
ed_from_gp_rate         = 1 / 2 # Lag from GP visit to ED visit
hosp_admit_direct_rate  = 1 / 4 # Lag from symptom onset to hospital admission (direct). Docherty et al. (2020); Knock et al. (2021): mean 4 days
hosp_admit_from_gp_rate = 1 / 1 # Lag from GP visit to hospital admission
```

### 12.2 Hospital general ward rates
Source: Knock et al. 2021, Table S2

```julia
hosp_recovery_rate            = 1 / 10.7  # General ward → recovery: 10.7d (95% CI: 0.3–39.4). Erlang(k=1, γ=0.09)
hosp_short_stay_recovery_rate = 1 / 0.49  # Short stay (<24h) mean duration. Estimated after splitting Erlang distribution at 24hrs.
hosp_long_stay_recovery_rate  = 1 / 11.70 # Long stay (≥24h) mean duration. Estimated after splitting Erlang distribution at 24hrs.
hosp_death_rate               = 1 / 10.3  # General ward → death: 10.3d (95% CI: 1.3–28.8). Erlang(k=2, γ=0.19)
```

### 12.3 ICU rates 
Source: Knock et al. 2021, Table S2

```julia
triage_icu_rate                          = 1 / 2.5  # Triage to ICU: 2.5d (95% CI: 0.1–9.2). Erlang(k=1, γ=0.4)
icu_to_death_rate                        = 1 / 11.8 # ICU → death: 11.8d (95% CI: 1.4–32.9). Erlang(k=2, γ=0.17)
icu_to_stepdown_leading_to_recovery_rate = 1 / 15.6 # ICU → stepdown (recovery): 15.6d (95% CI: 0.4–57.6). Erlang(k=1, γ=0.06)
stepdown_to_recovery_after_icu_rate      = 1 / 12.2 # Stepdown → recovery: 12.2d (95% CI: 1.5–34.0). Erlang(k=2, γ=0.16)
icu_to_stepdown_leading_to_death_rate    = 1 / 7.0  # ICU → stepdown (death): 7.0d (95% CI: 0.2–25.7). Erlang(k=1, γ=0.14)
stepdown_to_death_after_icu_rate         = 1 / 8.1  # Stepdown → death: 8.1d (95% CI: 0.2–29.7). Erlang(k=1, γ=0.12)
```

### 12.4 Discharge time limits

```julia
tdischarge_ed_upper_limit              = 0.5 # days — max time in ED before discharge (not admitted)
tdischarge_hosp_short_stay_upper_limit = 1.0 # days — max time for short-stay hospital admission
```

---

## 13. Sampling Parameters

### 13.1 ICU sampling

```julia
icu_sample_type               = "regional" # "regional" or "fixed". If "fixed" then {p_sampled_icu} function will be used, which doesn't take into account test sensitivity or practical sampling proportion, which "regional" option does by using {sample_icu_cases} function.
icu_site_stage                = "current"
p_sampled_icu                 = 0.15       # ICU sampling proportion (if icu_sample_type = "fixed")
sample_target_prob_icu        = 0.90       # Practical sampling proportion at ICU site
n_icu_samples_per_week        = 300        # Target number of ICU samples per week
icu_ari_admissions            = 1440       # Weekly ICU ARI admissions [793 summer, 1440 winter] estimated using unpublished UKHSA data
from the NHS Digital Secondary Users Survey
icu_ari_admissions_adult_p    = 0.76       # Proportion of ICU ARI admissions that are adults (≥16y) estimated using unpublished UKHSA data
from the NHS Digital Secondary Users Survey
icu_ari_admissions_child_p    = 0.24       # Proportion that are children (<16y) estimated using unpublished UKHSA data
from the NHS Digital Secondary Users Survey
turnaroundtime_icu            = [2, 4]     # [lower, upper] days, Uniform distribution
icu_swab_lag_max              = 1          # Max days between ICU admission and swab
icu_only_sample_before_death  = true       # Constrain swab time to before death to time of death if earlier than swab time
icu_nhs_trust_sampling_sites_file = "data/nhs_trust_site_sample_targets.csv"
```

### 13.2 Primary care sampling (based on the Oxford-RCGP RSC network)

```julia
turnaroundtime_rcgp = [2, 4]       # [lower, upper] days — swab to result. Source: UKHSA (2025a), Data quality report: national flu/COVID-19 surveillance report
gp_practices_total  = 6199         # Total GP practices in England. Source: BMA (2025) / NHS Digital (2025)
gp_practices_swab   = 300          # GP practices taking virology surveillance swabs
gp_swabs_mg         = 300          # Assumed number of swabs metagenomic sequenced
pop_eng             = 5.7106398e7  # Population of England. Source: UK ONS (2022) mid-year estimate
gp_ari_consults     = 327          # ARI consultations per 100k/week [180 summer, 327 winter 2024/25]. Estimated from Oxford-RCGP Research & Sureveillance Centre (RSC) (2025)
gp_ari_swabs        = 747          # Swabs from suspected ARI per week [319 summer, 747 winter 2024/25]. Estimated from Oxford-RCGP Research & Sureveillance Centre (RSC) (2025)
```

### 13.3 Secondary care sampling (based on the HARISS network)

```julia
turnaroundtime_hariss                = [2, 4]    # [lower, upper] days, Uniform distribution
hariss_courier_to_analysis           = 1.0       # Days from courier collection to analysis start
n_hosp_samples_per_week              = 300       # Total hospital samples per week
sample_allocation                    = "equal"   # "equal" or "weighted"
sample_proportion_adult              = "free"    # "free" or numeric proportion (e.g. 0.75)
weight_samples_by                    = "ae_mean" # or "catchment_pop". NHS Trust proportion of A&E attendances or NHS Trust catchment area population
phl_collection_dow                   = [2, 5]   # Day(s) of week for PHL courier collection (Mon=2, Thu=5)
phl_collection_time                  = 0.5      # Time of day for collection (0.5 = midday)
hosp_to_phl_cutoff_time_relative     = 1        # Days: Assume cutoff for swab to reach Public Health Laboratory before courier collection
swab_time_mode                       = 0.25     # Days: Assume swabbing peaks at 6h (0.25d) after admission
swab_proportion_at_48h               = 0.9      # Assume 90% of swabs taken within 48h of admission
proportion_hosp_swabbed              = 0.9      # Assume x% of ARI attendances are swabbed
hariss_only_sample_before_death      = true     # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
hariss_nhs_trust_sampling_sites_file = "data/hariss_nhs_trust_sampling_sites.csv"
```

### 13.4 Hospital background ARI admissions (for surveillance denominator)
Estimated from unpublished UKHSA ED Syndromic Surveillance System data: winter (Dec 2023, Jan 2024, Feb 2025) and summer (Jun, Jul, Aug 2025)
```julia
# ED ARI admissions
hosp_ari_admissions         = 6088 # Weekly ARI admissions [3452 summer]
hosp_ari_admissions_adult_p = 0.52 # Proportion adults (≥16y) [0.58 summer]
hosp_ari_admissions_child_p = 0.48 # Proportion children (<16y) [0.42 summer]

# ED ARI destination proportions (adults)
ed_ari_destinations_adult_p_discharged  = 0.628
ed_ari_destinations_adult_p_short_stay  = 0.030
ed_ari_destinations_adult_p_longer_stay = 0.342

# ED ARI destination proportions (children)
ed_ari_destinations_child_p_discharged  = 0.861
ed_ari_destinations_child_p_short_stay  = 0.014
ed_ari_destinations_child_p_longer_stay = 0.125
```

---

## 14. Metagenomic Test Sensitivity

Sensitivity of metagenomic testing by pathogen type. Source: Alcolea-Medina et al. (2025).

```julia
pathogen_type           = "virus" # Default pathogen type
sensitivity_mg_virus    = 0.89    # Sensitivity for viral pathogens
sensitivity_mg_bacteria = 0.97    # Sensitivity for bacterial pathogens
sensitivity_mg_fungi    = 0.89    # Sensitivity for fungal pathogens
```

---

## References

Alcolea-Medina et al. (2025), Rapid pan-microbial metagenomics for pathogen
detection and personalised therapy in the intensive care unit: a single-centre
prospective observational study, *Lancet Microbe*, 6(10). DOI:10.1016/j.lanmic.2025.101174

Bi et al. (2020), Epidemiology and transmission of COVID-19 in 391 cases and
1286 of their close contacts in Shenzhen, China: a retrospective cohort study, *Lancet Infectious Diseases*,
20(8). DOI: 10.1016/S1473-3099(20)30287-5

British Medial Association (2025). Pressures in general practice. Accessed: 2025-09-02. Available
from: https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pre
ssures/pressures-in-general-practice-data-analysis.

Charalampous et al. (2024), Routine Metagenomics Service for ICU Patients with
Respiratory Infection, *Am J Respir Crit Care Med*, 209(2), pp 164–174.
DOI: 10.1164/rccm.202305-0901OC

Chen et al. (2022), Inferring time-varying generation time, serial interval, 
and incubation period distributions for COVID-19, *Nature Communications*, 13, 7727.
DOI: 10.1038/s41467-022-35496-8

Danon et al. (2013), Social encounter networks: characterizing Great Britain,
*Proceedings of the Royal Society B: Biological Sciences*, 280(1765).
DOI: 10.1098/rspb.2013.1037

Delignette-Muller & Dutang (2015), fitdistrplus: An R Package for Fitting
Distributions, *Journal of Statistical Software*, 64(4), 1–34.
DOI: 10.18637/jss.v064.i04

Docherty et al. (2020), Features of 20,133 UK patients in hospital with
COVID-19 using the ISARIC WHO Clinical Characterisation Protocol, *BMJ*,
369:m1985. DOI: 10.1136/bmj.m1985

du Plessis et al. (2020), Establishment and lineage dynamics of the SARS-CoV-2
epidemic in the UK, *Science*, 371(6530), pp 708–712.
DOI: 10.1126/science.abf2946

Galmiche et al. (2023), SARS-CoV-2 incubation period across variants of
concern, *Lancet Microbe*, 4(6), e409–e417.
DOI: 10.1016/S2666-5247(23)00005-8

Hart et al. (2022), Inference of the SARS-CoV-2 generation time using UK household data,
*eLife*, 11. DOI: 10.7554/eLife.70767

Knock et al. (2021), Key epidemiological drivers and impact of interventions in
the 2020 SARS-CoV-2 epidemic in England, *Sci. Transl. Med.*, 13, eabg4262.
DOI: 10.1126/scitranslmed.abg4262

Lau et al. (2021), Joint Estimation of Generation Time and Incubation Period for Coronavirus Disease 2019, 
*The Journal of Infectious Diseases*, 224(10). DOI: 10.1093/infdis/jiab424

Mossong et al. (2008), Social Contacts and Mixing Patterns Relevant to the
Spread of Infectious Diseases, *PLoS Medicine*, 5(3), e74.
DOI: 10.1371/journal.pmed.0050074.
Survey data: https://zenodo.org/records/3874557

NHS England (2025). NHS Digital General Practice Workforce Statistics. Accessed: 2025-09-02.
Available from: https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services.

Oxford-RCGP Research & Sureveillance Centre (RSC) (2025). Virology Dashboard. Accessed: 2025-08-25. Available from: https://orchid.phc.ox.ac.uk/surveillance/dashboards-and-observatories-portal/virology-dashboard.

Reddy et al (2026), A multi-scale model to evaluate airport wastewater surveillance and ICU genomic monitoring for pandemic preparedness, *medRxiv*, DOI: 10.64898/2026.02.27.26347250

Saigal et al. (2025), Predictors of specialist care referrals (SCR) following emergency department review or hospital admission in adults with previous acute COVID-19: a prospective UK cohort study, *BMC Emergency Medicine*, 25(1), DOI: 10.1186/s12873-024-01164-x

St-Onge et al. (2025), Pandemic monitoring with global aircraft-based wastewater surveillance networks, *Nature Medicine*, 2025;31:788–779, DOI: 10.1101/2024.08.02.24311418

UK Health Security Agency (2025a), Data quality report: national flu and COVID-19 surveillance report, https://www.gov.uk/government/publications/sources-of-surveillance-data-for-influenza-covid-19-and-other-respiratory-viruses/data-quality-report-national-flu-and-covid-19-surveillance-report#primary-care-surveillance, [Accessed: 2025-05-27]

UK Health Security Agency (2025b), Flu Survey 2024/25 - Influenza in the UK, annual epidemiological report: winter 2024 to 2025, https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025, [Accessed: 2025-09-04]

UK Office for National Statistics (2022). Estimates of the population for the UK, England, Wales,
Scotland, and Northern Ireland - Mid-2022 Edition. Accessed: 2024-11-06. Available from:
https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/p
opulationestimates/datasets/populationestimatesforukenglandandwalesscotlanda
ndnorthernireland.

Verity et al. (2020), Estimates of the severity of coronavirus disease 2019: a model-based analysis,
*The Lancet Infectious Diseases*, 20(6).
DOI: 10.1016/S1473-3099(20)30243-7

Volz et al. (2021), Evaluating the Effects of SARS-CoV-2 Spike Mutation D614G
on Transmissibility and Pathogenicity, *Cell*, 184(1), 64–75.
DOI: 10.1016/J.CELL.2020.11.020

WHO (2020). Report of the WHO-China Joint Mission on Coronavirus Disease 2019
(COVID-19). Feb 28, 2020.

Xu et al. (2023), Assessing changes in incubation period, serial interval, and generation time of SARS-CoV-2 variants of concern: a systematic review and meta-analysis, *BMC Medicine*, 21:374.
DOI: 10.1186/s12916-023-03070-8

Zhao et al. (2021), Estimating the generation interval and inferring the latent period of COVID-19 from the contact tracing data, *Epidemics*, 36. DOI:10.1016/j.epidem.2021.100482.

---
