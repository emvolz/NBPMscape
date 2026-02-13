### Load severity probabilities
# Extracted from Knock et al (2021), Key Epidemiological Drivers and Impact of Interventions in the 2020 SARS-CoV-2 Epidemic in England’.
# Science Translational Medicine, 16(602), eabg4262
symptomatic_prob_by_age = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_symptomatic_prob_by_age.rds" ) )
symptomatic_prob_by_age = filter(row -> !(row.age_single_year in ["Care_home_workers", "Care_home_residents"])
                                      , symptomatic_prob_by_age)
const SYMPTOMATIC_PROB_BY_AGE = symptomatic_prob_by_age[:,["age_single_year","symptomatic_prob"]]
ihr_ifr_by_age = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_ihr_ifr_by_age.rds" ) )
ihr_ifr_by_age = filter(row -> !(row.age_single_year in ["Care_home_residents"])
                      , ihr_ifr_by_age)
ihr_ifr_by_age[:,(3:8)] = ihr_ifr_by_age[:,(3:8)] ./ 100 # Convert from percentage to proportion
const IHR_BY_AGE = ihr_ifr_by_age[:,["age_single_year","IHR"]]
const IFR_BY_AGE = ihr_ifr_by_age[:,["age_single_year","IFR"]]

#const IFR_BY_AGE = parse.(Float64, ihr_ifr_by_age.age_single_year)","IFR"])
#care_pathway_prob_by_age
const CARE_PATHWAY_PROB_BY_AGE = load( joinpath( @__DIR__, "..", "data/severity_care_probabilities", "knock_care_pathway_prob_by_age.rds" ) )
