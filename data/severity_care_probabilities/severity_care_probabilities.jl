### Load severity probabilities

# Load, manipulate and save COVID-19 parameters from Knock et al (2021) in format with the 
# minimal data required by parameter configuration method.
# Extracted from Knock et al (2021), Key Epidemiological Drivers and Impact of Interventions in the 2020 SARS-CoV-2 Epidemic in England’.
# Science Translational Medicine, 16(602), eabg4262
# Symptomatic probability by age
symptomatic_prob_by_age = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_symptomatic_prob_by_age.rds" ) )
symptomatic_prob_by_age = filter(row -> !(row.age_single_year in ["Care_home_workers", "Care_home_residents"])
                                      , symptomatic_prob_by_age)
symptomatic_prob_by_age = symptomatic_prob_by_age[:,["age_single_year","symptomatic_prob"]]
# Infection hospitalisation rates (IHR) and infection fatality rates (IFR) by age
ihr_ifr_by_age = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_ihr_ifr_by_age.rds" ) )
ihr_ifr_by_age = filter(row -> !(row.age_single_year in ["Care_home_residents"])
                      , ihr_ifr_by_age)
ihr_ifr_by_age[:,(3:8)] = ihr_ifr_by_age[:,(3:8)] ./ 100 # Convert from percentage to proportion
ihr_ifr_by_age = ihr_ifr_by_age[:,["age_single_year","IHR","IFR"]]
# Merge symptomatic probability, IHR and IFR by age
symptomatic_IHR_IFR_by_age = outerjoin( symptomatic_prob_by_age, ihr_ifr_by_age, on = :age_single_year )
# Save to .csv
CSV.write( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "covid_knock_symptomatic_ihr_ifr_by_age.csv" ) 
          , symptomatic_IHR_IFR_by_age )

# Care pathway probabilities by age
care_pathway_prob_by_age = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_care_pathway_prob_by_age.rds" ) )
care_pathway_prob_by_age = care_pathway_prob_by_age[:, ["age_single_year","p_hosp_sympt","p_ICU_hosp","p_death_ICU","p_death_hosp_D","p_death_stepdown"]]
CSV.write( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "covid_knock_care_pathway_prob_by_age.csv" ) 
          , care_pathway_prob_by_age )