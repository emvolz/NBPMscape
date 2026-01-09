Plots.default(;lw = 4)

# Include configuration functions
include("config.jl")

const MAXDURATION = 180.0

mutable struct Infection
	pid::String # infection id 
	dpid::Union{String,Missing}
	H::DataFrame # transmission linelist 
	R::Int64 # number transmissions 
	sol::Union{ODESolution,Nothing}
	tinf::Float64 # time infected 
	tgp::Float64 # time visited GP
	ted::Float64 # time visited Emergency Department (ED)
	thospital::Float64 # time admitted to hospital general ward
	ticu::Float64 # time admitted to ICU
	tstepdown::Float64 # time stepped down from ICU to general ward
	tdischarge::Float64 # time discharged from hospital (general ward, ICU, or stepdown)
	trecovered::Float64 # time recovered - either community, discharge from hospital, ICU, or stepdown ward
	tdeceased::Float64 # time of death - either in hospital, ICU, or stepdown ward
	contacttype::Symbol # cause of infection
	degree::Tuple{Int64,Int64,Float64}
	iscommuter::Bool
	initialregion::String
	homeregion::String
	commuteregion::String
	initialdayofweek::Int # day 1-7 when infection ocurred 
	severity::Symbol #severity::NamedTuple{(:severity, :fatal), Tuple{Symbol, Bool}} # severity of infection and whether it will be fatal
	fatal::Bool # Will the infection be fatal {true,false}. Only possible to be true for severity = verysevere or severe
	generation::Int64
	isdeceased::Bool # Whether individual has died (default = FALSE). Assume all individuals that die have been admitted to ICU. 
	donor_age::Union{Int8,Missing} # Age of infector
	infectee_age::Union{Int8,Missing} # using single year age so can incorporate different age groupings for different data inputs, e.g. contact number distributions may have different age groups to the number of ICU admissions
	importedinfection::Bool # Record whether individual was infected when arrived in UK - if they are admitted to ICU, this may prompt further investigation
end


const CARE  = (:undiagnosed, :GP, :ED, :admittedhospital, :admittedicu, :stepdown, :discharged, :deceased)
const SEVERITY = [:asymptomatic, :mild, :moderate_GP, :moderate_ED, :severe_hosp_short_stay, :severe_hosp_long_stay, :verysevere ] #[:asymptomatic, :mild, :moderate, :severe, :verysevere ] #const SEVERITY = [:mildorasymptomatic, :moderate, :severe, :verysevere ]
const STAGE = (:latent, :infectious, :recovered, :deceased)

const CONTACTTYPES = (:F, :G, :H)

## Function to convert age range strings (e.g. "0-18") to Range type (e.g. 0:18)
function parse_range(str::String)
	parts = parse.(Int, split(str, ':'))
	return parts[1]:parts[2]
end							 

## Define contact number distributions, stratified by age of individual and contact type
contact_age_groups = replace.( sort( unique( CONTACT_DISTRIBUTIONS.age_group ) , rev = false )
							 , "_" => ":", "plus" => ":100")  # amend the format of the age group descriptions (i.e. separator and "75plus" -> 75:100)
contact_age_groups_range = parse_range.( contact_age_groups )  # Convert from String to Range type

# Build dataframe containing contact number distributions by age group
contact_rate_dist_par_age_groups = DataFrame(
    age_group    = contact_age_groups_range
	# f is household size distribution 
	,fnegbinomr   = sort( filter( row -> row.contact_setting =="home", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    ,fnegbinomp   = sort( filter( row -> row.contact_setting =="home", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	# g is Work / School contact number distribution
	# POLYMOD approximate work or school contact number distribution
	
    ,gnegbinomr   = sort( filter( row -> row.contact_setting =="work_school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    ,gnegbinomp   = sort( filter( row -> row.contact_setting =="work_school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
    
	# Other contact number distribution
	# POLYMOD approximate contact number distributions for combined settings: transport, leisure, otherplace
    ,oorateshape  = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    ,ooratescale  = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
    # excess rate distribution (also gamma) oorateshape1 = oorateshape + 1 and ooratescale1 = ooratescale
	,oorateshape1 = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape .+ 1
	,ooratescale1 = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
)

# Create expanded contact number distribution parameter vectors
max_age = maximum(last.(contact_rate_dist_par_age_groups.age_group))
contact_dist_params = [:fnegbinomr, :fnegbinomp
					 , :gnegbinomr, :gnegbinomp
					 #, :trateshape, :tratescale
					 , :oorateshape, :ooratescale, :oorateshape1, :ooratescale1]

# Generate vectors of contact number distribution parameters for each single year age up to maximum age
# and store in dictionary
contact_rate_dist_par_all_ages = Dict( param => Vector{Float64}(undef, max_age + 1) for param in contact_dist_params)
# Assign distribution parameter values to single age elements
for row in eachrow( contact_rate_dist_par_age_groups )
    age_range = row.age_group
    for param in contact_dist_params
		# Note that the first element of each vector is for age = 0 years
        contact_rate_dist_par_all_ages[param][age_range .+ 1] .= row[param] 
    end
end

### Define assortativity for contact ages

# Function to convert contact matrix by age group to single year age
function cont_matrix_age_group_to_single_yr(; contact_matrix_by_age_group = :nothing
											, single_year_ages = :nothing)

	# Initialise single year age contact matrix
	contact_matrix_age = NamedArray( zeros( maximum(single_year_ages)+1, maximum(single_year_ages)+1)
									,( single_year_ages, single_year_ages)
									)
	# Fill elements of single year age contact matrix
	for (row_idx, age_range_row) in enumerate(names(contact_matrix_by_age_group, 1)),
		(col_idx, age_range_column) in enumerate(names(contact_matrix_by_age_group, 2))

		value = contact_matrix_by_age_group[row_idx, col_idx]
		row_range = parse_range(age_range_row)
		col_range = parse_range(age_range_column)
		
		for r in row_range, c in col_range
			contact_matrix_age[r+1, c+1] = value
		end
	end
	return( contact_matrix_age )
end

single_year_ages = 0:100

const CONTACT_MATRIX_HOME_SINGLE_YEAR = cont_matrix_age_group_to_single_yr( contact_matrix_by_age_group = CONTACT_MATRIX_HOME
																			, single_year_ages = single_year_ages )
const CONTACT_MATRIX_SCHOOL_WORK_SINGLE_YEAR = cont_matrix_age_group_to_single_yr( contact_matrix_by_age_group = CONTACT_MATRIX_SCHOOL_WORK
																			, single_year_ages = single_year_ages )
const CONTACT_MATRIX_OTHER_SINGLE_YEAR = cont_matrix_age_group_to_single_yr( contact_matrix_by_age_group = CONTACT_MATRIX_OTHER
																			, single_year_ages = single_year_ages )



"""
Function    initialize_parameters(config_file::String="")

Decription	Initialize the global parameter structure P. Location of .yaml file containing the required
			parameters should be input as a string. If no argument is input then default parameters are used.

Argument	config_file::String		Path and filename for .yaml file containing model parameters.

Example 	# Note that ending the line with a semi-colon is important to stop the full set of parameters 
			# being printed to the screen
			initialize_parameters(); # Returns default parameters
			or
			initialize_parameters("config/outbreak_params_covid19_like.yaml"); # Updates parameters from file
			config_file = "config/outbreak_params_covid19_like.yaml"
"""
function initialize_parameters(config_file::String="")
	global P;
	# Reset P so if there is a failure in loading the new parameters there
	# is no mistake with an old parameter set remaining
	P = nothing 
    # Create default parameters
	# There are two sets: data-dependent and configurable. 
	# The latter can be replaced with values defined in the config_file while the former are not.
    #P = create_default_parameters()
	default_params = create_default_parameters(); # split(default_params)
	# Also define the default configurable parameters so can check 
	# them against the values coming from the config file
	default_configurable_params = default_params.default_configurable_params; # typeof(default_configurable_params); collect(propertynames(default_configurable_params))
    # Apply configuration if provided
    if !isempty(config_file)
        #try
            config_data = load_config(config_file);
            warnings, errors = validate_config(config_data);
            # Process warnings
			for warning in warnings
                @warn warning
            end
			# If there are validation errors, throw immediately
            if !isempty(errors)
                for error in errors
                    @error error
                end
                throw(ArgumentError("Configuration validation failed"))
            end
			# Only proceed if validation is passed
            P = update_configurable_parameters(default_params.default_P, config_data, default_configurable_params); #update_configurable_parameters(P, config_data, default_configurable_params);
            @info "Successfully loaded configuration from: $config_file"
			# Print changes made to configurable parameters during parameter initialization
			print_changes(P, config_data, default_configurable_params)
        #catch e
		#	# Only catch specific non-validation errors for graceful fallback
		#	if isa(e, ArgumentError) && occursin("validation failed", e.msg)
        #        # Re-throw validation errors - don't handle 'gracefully'
        #        rethrow(e)
        #    elseif isa(e, SystemError) || isa(e, Base.IOError)
        #        # File not found errors - handle 'gracefully' and use default parameters
        #        @error "Failed to load configuration from $config_file: $e"
        #        @info "Using default parameters"
        #    elseif isa(e, YAML.ParserError) || isa(e, YAML.ScannerError)
        #        # YAML parsing errors - handle 'gracefully' and use default parameters
        #        @error "Failed to parse YAML configuration file $config_file: $e"
        #        @info "Using default parameters"
        #    else
        #        # For validation errors and other critical errors, re-throw
        #        rethrow(e)
        #    end
        #end
    else
        @info "Using default parameters (no config file specified)"
		P = default_params.default_P;
    end
	# Use some parameter values to fill dataframes
	P = convert_params_to_dfs(P);
    return P;
end

function create_default_parameters()
    # Define parameters that depend on data objects
	data_dependent_params = (
		# Contact rate distributions stratified by age of individual and contact type
		# Based on POLYMOD UK - Mossong et al (2008), PLOS Medicine, 5(3)
		# Note that element 1 is for age = 0 years
		  fnegbinomr   = contact_rate_dist_par_all_ages[:fnegbinomr]
		, fnegbinomp   = contact_rate_dist_par_all_ages[:fnegbinomp]
		, gnegbinomr   = contact_rate_dist_par_all_ages[:gnegbinomr]
		, gnegbinomp   = contact_rate_dist_par_all_ages[:gnegbinomp]
		, oorateshape  = contact_rate_dist_par_all_ages[:oorateshape]
		, ooratescale  = contact_rate_dist_par_all_ages[:ooratescale]
		, oorateshape1 = contact_rate_dist_par_all_ages[:oorateshape1]
		, ooratescale1 = contact_rate_dist_par_all_ages[:ooratescale1]

		# Contact rate matrix by age
		# Based on POLYMOD UK - Mossong et al (2008), PLOS Medicine, 5(3)
		, f_contact_matrix_age = CONTACT_MATRIX_HOME_SINGLE_YEAR
		, g_contact_matrix_age = CONTACT_MATRIX_SCHOOL_WORK_SINGLE_YEAR
		, o_contact_matrix_age = CONTACT_MATRIX_OTHER_SINGLE_YEAR

		### Infection severity disaggregated by age
		
		## Single values calculated using disaggregated probabilities (see below) and
		## weighted by population age distribution
		#, propverysevere = 0.008 # Proportion of infected that are admitted to ICU. Based on age group disaggregated probabilities in Knock et al (2021) weighted by ONS England population data for MYE 2022 by age group. Also supported by value reported in Thygesen et al (2022) Lancet Digital Health, 6.4% admitted to hospital, of which 10.6% admitted to ICU = 0.7%
		#, propsevere = 0.03 # Proportion of infected that are admitted to hospital but not ICU. Based on age group disaggregated probabilities in Knock et al (2021) weighted by ONS England population data for MYE 2022 by age group. Also supported by value reported in Thygesen et al (2022) Lancet Digital Health, 6.4% admitted to hospital, of which 10.6% admitted to ICU = 0.7%
		#, propmoderate = 0.053 # Analysis of FluSurvey data for 2024/25 shows 11% (weekly mean) of symptomatic ILI consult GP in person or via phone. Symptomatic proportion is 49.3%.
		#, prop_asymptomatic = 0.507 # = 1 - weighted mean of symptomatic probability in Knock et al (2021) (weighted by MYE 2022 ONS population of England by single year age)
		#, propmildorasymptomatic = 0.507 + (1 - 0.008 - 0.03 - 0.053 - 0.507) # prop_asymptomatic + (1-prop_asymptomatic-propmoderate-propsevere-propverysevere)  # asymptomatic + mild (a balancing number), i.e. no healthcare and no sampling
		# CHECK: propverysevere + propsevere + propmoderate + propmildorasymptomatic == 1
		
		## Infection severity probabilities disaggregated by age
		## Sourced from Knock et al (2021), Science Translational Medicine, 13(602)
		, symptomatic_prob_by_age = SYMPTOMATIC_PROB_BY_AGE
		, ihr_by_age = IHR_BY_AGE[:,"IHR"]
		, symptomatic_ihr_by_age = parse.(Float64, CARE_PATHWAY_PROB_BY_AGE[:,"p_hosp_sympt"])
		, icu_by_age = parse.(Float64, CARE_PATHWAY_PROB_BY_AGE[:,"p_ICU_hosp"]) # Prob of admission to ICU if already admitted to hospital

		## Probability of death by age and care stage
		## Sourced from Knock et al (2021), Science Translational Medicine, 13(602)
		, ifr_by_age = IFR_BY_AGE[:,"IFR"] # Infection fatality ratio
		, p_death_icu = parse.(Float64, CARE_PATHWAY_PROB_BY_AGE[:,"p_death_ICU"])
		, p_death_hosp = parse.(Float64, CARE_PATHWAY_PROB_BY_AGE[:,"p_death_hosp_D"])
		, p_death_stepdown = parse.(Float64, CARE_PATHWAY_PROB_BY_AGE[:,"p_death_stepdown"])
		#plot(ifr_by_age); plot!(p_death_icu); plot!(p_death_hosp);plot!(p_death_stepdown)

		# Hospital parameters
        , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
        #, nhs_trust_ae_12m = AE_12M
	);

    
    # Define configurable parameters
	configurable_params = ( 
		# Relative contact rates
		# Normalised to 1 from the 10.76 contact hours in Table S2 of 
		# Danon et al (2013), Proceedings of the Royal Society B: Biological Sciences, 280(1765)
		fcont = 10.76 / 10.76   # for flinks (household) 
		, gcont = 6.71 / 10.76  # for glinks (work/school)
		, oocont = 8.09 / 10.76 # for other

		# Scale adjustment for contact rates by day of the week (Sun-Sat)
		# Sourced from POLYMOD - Mossong et al (2008), PLOS Medicine, 5(3)
		#           Sun      , Mon      , Tue      , Wed      , Thur     , Fri      , Sat 
		, dowcont = (0.1043502, 0.1402675, 0.1735913, 0.1437642, 0.1596205, 0.1445298, 0.1338766) 
		# Sourced from 'Social Contact Survey' 
		# Danon et al (2013), Proceedings of the Royal Society B: Biological Sciences, 280(1765)
		# Table S1 mean degree
		#           Sun      , Mon      , Tue      , Wed      , Thur     , Fri      , Sat 
		#, dowcont = ([28.00,25.39,30.53,27.21,26.82,28.80,21.76]/sum([28.00,25.39,30.53,27.21,26.82,28.80,21.76])) 
		# [0.1485, 0.1347, 0.1620, 0.1443, 0.1423, 0.1528, 0.1154]
		# Table S1 contact hours
		#           Sun      , Mon      , Tue      , Wed      , Thur     , Fri      , Sat 
		#, dowcont = ([25.80,24.51,26.47,27.40,26.72,26.50,26.46]/sum([25.80,24.51,26.47,27.40,26.72,26.50,26.46])) 
		# [0.1403, 0.1333, 0.1439, 0.1490, 0.1453, 0.1441, 0.1439]
		
		, infectivity = 2.00 #1.25 # scales transmission rate # Updated to 2.00 after incorporation of age disaggregated parameters # Use infectivitytoR() to check R value for current parameters
		, infectivity_shape = 2.2 * 0.75 # # Manually calibrated to give a generation time (Tg) of 5-6 days given the infectious period distribution (parameters below - Verity et al (2020)). Hart et al (2022), eLife, 11, DOI: 10.7554/eLife.70767. Lau et al (2021), The Journal of Infectious Diseases, 224(10), DOI:10.1093/infdis/jiab424. Chen et al (2022), Nature, 13, 7727, DOI:10.1038/s41467-022-35496-8. Xu et al (2023), BMC Medicine, 21:374, DOI:10.1186/s12916-023-03070-8.  Bi et al (2020), Lancet Infectious Diseases, 20(8), DOI:10.1016/S1473-3099(20)30287-5. 
		, infectivity_scale = 2.5 * 0.75
		
		, latent_shape = 3.26 # Distribution parameters inferred in 'latent_period_estimate.R' from results reported in Zhao et al (2021), Epidemics, 36, 100482 - mean latent period of 3.3 days (95% CI 0.2, 7.9)
		, latent_scale = 0.979
		, infectious_shape = 8.16 # mean 24d from Verity et al (2020), The Lancet Infectious Diseases, 20(6)
		, infectious_scale = 3.03 

		, ρ_hosp = 0.250 #  transmission reduction, i.e. transmission rate is only 25% of normal
		, ρ_asymptomatic = 0.223 #  transmission reduction for asymptomatic individuals - see Knock et al (2021) Supplementary Information

		, frate = 0.0 # rate of gaining & losing flinks
		, grate = 1/30.0 # rate of gaining and losing

		# Severity probabilities with no age disaggregation
		# Split 'severe' infections, which is defined by hospital admission, into short and long stay.
		, prop_severe_hosp_short_stay = 0.40351 # Given admitted to hospital, will stay for <24h. Estimated using Knock et al (2021) and Saigal et al. (2025), BMC Emergency Medicine, 25:11. See 'HARISS_pathway_adjustment.jl' 
		, prop_severe_hosp_long_stay  = 0.59649 # Given admitted to hospital, will stay for >24h. Estimated using Knock et al (2021) and Saigal et al. (2025), BMC Emergency Medicine, 25:11. See 'HARISS_pathway_adjustment.jl' 
		# Note that mild, moderate_GP and moderate_ED sum to 1. They are symptomatic but the distinction is whether they do not seek healthcare, consult a GP (in person or via phone) or visit an Emergency Department (ED) respectively.
		, prop_moderate_ED = 0.12445 # should be approximately equivalent to 0.026 of all infections # Will visit hospital Emergency Department (ED) but will be discharged without admission. Estimated using Knock et al (2021) and Saigal et al. (2025), BMC Emergency Medicine, 25:11. See 'HARISS_pathway_adjustment.jl' 
		, prop_moderate_GP = 0.11790 # Will visit a GP but won't go to hospital. Estimated from FluSurvey data (people with ILI - not COVID specific), combining "Phoned GP" and "Visited GP" categories, for 2024/25 season (not full 12 month period) in Figure 11 of 'Influenza in the UK, annual epidemiological report: winter 2024 to 2025’, published on 22 May 2025 [Accessed on 4 Sep 2025 at https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025 
		, prop_mild     = 0.75765 # = 1 - prop_moderate_GP - prop_moderate_ED  # Symptomatic but won't visit a GP or go to hospital
		
		## Rates (= 1 / duration) from symptom onset to various care stages
		# Where a patient visits a GP before visiting hospital, the rates are such that the number of days is approximately the
		# same as going directly to hospital (either ED and discharge or hospital admission).
		# For those admitted to hospital, we do not add any time for attendance at ED before admission.
		, gp_only_rate = 1/5 # Typical lag of 3-7 days between ARI symptom onset (not pathogen specific) and GP visit (Data quality report: national flu and COVID-19 surveillance report (27 May 2025))
		, ed_direct_rate = 1/5 # Rate of visiting a hospital Emergency Department (ED) but not being admitted. Therefore assume to be the same as the GP rate.
		, gp_before_hosp_rate = 1/3 # Less than gp_only_rate to allow infected individual to reach hospital in approximately the same time as if went directly to hospital without visiting the GP.
		, ed_from_gp_rate = 1/2 # Rate of visiting a hospital Emergency Department (ED) but not being admitted. Therefore assume to be the same as the GP rate.
		, hosp_admit_direct_rate = 1/4 # Docherty et al (2020), BMJ, 369:m1985, doi: 10.1136/bmj.m1985 # Also in Knock et al (2021), Science Translational Medicine, 13(602) "mean time from symptom onset to admission to hospital" = 4 days
		, hosp_admit_from_gp_rate = 1/1 # combined with gp_before_hosp_rate should be approximately the same as hosp_admit_direct_rate
		
		## Rate (= 1 / duration) in different care stages for different pathways - source Knock et al (2021), Science Translational Medicine, 13(602), Table S2
		# General ward rates
		, hosp_recovery_rate = 1 / 10.7 # 'Hospitalised on general ward leading to recovery' 10.7 days (95% CI: 0.3-39.4). Erlang(k=1,gamma=0.09). Knock et al (2021). 
		, hosp_short_stay_recovery_rate = 1 / 0.49 # hosp_recovery_rate above split into short stay (<24h) and long stay (>24h). Mean durations of the two periods computed using erlang_truncated_means( k = 1, rate = 0.0935, lower_limit = 0.0, mid_limit = 1.0, upper_limit = Inf). Note that the rate has been adjusted as the rounded rate parameter does not give the mean value reported.
		, hosp_long_stay_recovery_rate  = 1 / 11.70 # hosp_recovery_rate above split into short stay (<24h) and long stay (>24h). Mean durations of the two periods computed using erlang_truncated_means( k = 1, rate = 0.0935, lower_limit = 0.0, mid_limit = 1.0, upper_limit = Inf). Note that the rate has been adjusted as the rounded rate parameter does not give the mean value reported.
		, hosp_death_rate    = 1 / 10.3 # 'Hospitalised on general ward leading to death' 10.3 days (95% CI: 1.3-28.8). Erlang(k=2,gamma=0.19). Knock et al (2021). 
		# ICU rates
		, triage_icu_rate    = 1 /  2.5 # 'Triage to ICU' 2.5 days (95% CI: 0.1-9.2). Erlang(k=1,gamma=0.4). Knock et al (2021). 
		, icu_to_death_rate  = 1 / 11.8 # 'Hospitalised in ICU, leading to death' 11.8 days (95% CI: 1.4-32.9). Erlang(k=2,gamma=0.17). Knock et al (2021). 
		, icu_to_stepdown_leading_to_recovery_rate = 1 / 15.6 # 'Hospitalised in ICU, leading to recovery' 15.6 days (95% CI: 0.4, 57.6). Erlang(k=1, gamma=0.06) . Knock et al (2021)
		# ICU stepdown rates
		, stepdown_to_recovery_after_icu_rate  = 1 / 12.2 # 'Stepdown recovery period after leaving ICU' 12.2 days (95% CI: 1.5-34.0). Erlang(k=2,gamma=0.16). Knock et al (2021). 
		
		, tdischarge_ed_upper_limit = 0.5 # days. People attending Emergency Department (ED) (and not admitted) will be discharged by this time.
		, tdischarge_hosp_short_stay_upper_limit = 1.0 # days. People admitted to hospital for a short stay will be discharged by this time.

		, commuterate = 2.0

		, importrate = .5 # if using constant rate 
		, nimports = 1000 # du plessis 2020
		# t distribution
		, import_t_df = 2.48 # personal analysis of lineages studied in Volz et al. Cell 2020 
		, import_t_s = 8.75
	#      df           m           s    
	#   2.482146   -1.074888    8.749736 
	# ( 1.226155) ( 1.608095) ( 1.978287)

		, μ = 0.001 # mean clock rate -- additive relaxed clock
		, ω = 0.5 # variance inflation

		# Sampling parameters
		# ICU
		, icu_sample_type = "regional" # "regional" or "fixed". If "fixed" then p_sampled_icu will be used, note that this doesn't take into account test sensitivity or practical sampling proportion, which "regional" does by using {sample_icu_cases} function
		, icu_site_stage = "current" #
		, sample_icu_cases_version = "number" # "proportion". Is sampling quantity based on a proportion of all cases or a fixed number of samples?
		, p_sampled_icu = 0.15 # ICU sampling proportion
		, sample_target_prob_icu = 0.90 # To account for some patients not being tested due to logistics/practicalities even though ideally 100% are tested at a particular ICU site
		, n_icu_samples_per_week = 300 # Number of samples to be taken from ICU per week
		, icu_ari_admissions = 1440 # [793, 1440] # Weekly ICU admission numbers [summer,winter]
		, icu_ari_admissions_adult_p = 0.76 # Proportion of ICU ARI admissions that are adults (16y and over)
		, icu_ari_admissions_child_p = 0.24 # Proportion of ICU ARI admissions that are children (<16y)
		, turnaroundtime_icu = [2,4] # upper and lower limits, in days, of time taken to process sample and report results /declare detection. tsample drawn from a Uniform(lower, upper) distribution in sampling functions
		, icu_swab_lag_max = 1 # days. Upper limit on the time between admission to ICU and a swab being taken. Simulated time = Uniform(ticu, ticu + icu_swab_lag_max)
		, icu_nhs_trust_sampling_sites_file = "data/nhs_trust_site_sample_targets.csv" #CSV.read( "data/nhs_trust_site_sample_targets.csv" , DataFrame )
		, icu_only_sample_before_death = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
		# Primary care (i.e. GPs via RCGP)
		, turnaroundtime_rcgp = [2,4] # upper and lower limits, in days, of time between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing. tsample drawn from a Uniform(lower, upper) distribution in sampling functions.
		, gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
		, gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
		, gp_swabs_mg = 300 # [319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
        , pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
        , gp_ari_consults = 327 # 180, 327 # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
        , gp_ari_swabs = 747 #[319, 747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
        # Secondary care (i.e. hospitals via HARISS network)
		, turnaroundtime_hariss = [2,4] # upper and lower limits, in days, of time taken to process sample and report results /declare detection. tsample drawn from a Uniform(lower, upper) distribution in sampling functions
		, initial_dow = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
		, hariss_courier_to_analysis = 1.0 # Time between courier collection from PHL and beginning of analysis
        , n_hosp_samples_per_week = 300 # Total number of hospital samples to be taken per week
        , sample_allocation = "equal" # "equal" or "weighted"
        , sample_proportion_adult = "free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target between adults and children. "free" indicates that no split is specified
        , hariss_nhs_trust_sampling_sites_file = "data/hariss_nhs_trust_sampling_sites.csv"#CSV.read("data/hariss_nhs_trust_sampling_sites.csv", DataFrame) # List of NHS Trusts in HARISS sampling network 
        , weight_samples_by = "ae_mean" # or "catchment_pop". NHS Trust proportion of A&E attendances or NHS Trust catchment area population
        , phl_collection_dow = [2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
        , swab_time_mode = 0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
        , swab_proportion_at_48h = 0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
        , proportion_hosp_swabbed = 0.9 # Assume X% of ARI attendances are swabbed
        , hariss_only_sample_before_death = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
        # Hospital parameters
        # Seasonal values
        # Winter
        , hosp_ari_admissions = Int64( round( 79148 / ((31+31+29)/7), digits = 0 ) ) # for winter and Int64(45360 / ((30+31+31)/7)) for summer. Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated) - using Dec 2023, Jan 2024, Feb 2025 data for winter and Jun, Jul, Aug 2024 for summer
        , hosp_ari_admissions_adult_p = 0.52 # for winter and 0.58 for summer.Proportion of ED ARI admissions that are adults (16y and over)
        , hosp_ari_admissions_child_p = 0.48 # for winter and 0.42 for summer. Proportion of ED ARI admissions that are children (<16y)
        #, ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
        #                                                           , proportion_of_attendances = [0.628,0.030,0.342])
        # , ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
        #                                                           , proportion_of_attendances = [0.861,0.014,0.125])
        , ed_ari_destinations_adult_p_discharged  = 0.628
		, ed_ari_destinations_adult_p_short_stay  = 0.030
		, ed_ari_destinations_adult_p_longer_stay = 0.342
        , ed_ari_destinations_child_p_discharged  = 0.861
		, ed_ari_destinations_child_p_short_stay  = 0.014
		, ed_ari_destinations_child_p_longer_stay = 0.125
		
		, pathogen_type = "virus"		

		# Sensitivity of metagenomic testing 
		# see Alcolea-Medina et al (2025), "Rapid pan-microbial metagenomics for pathogen detection and personalised therapy in the intensive care unit: a single-centre prospective observational study", Lancet Microbe, 
		, sensitivity_mg_virus = 0.89
		, sensitivity_mg_bacteria = 0.97
		, sensitivity_mg_fungi = 0.89

		## Experimental features yet to be incorporated into model
		## Assumptions regarding current non-metagenomic surveillance
		# Proportion of ICU admissions with history of international travel that would be 
		# further investigated (with non-metagenomic surveillance)
		, non_mg_inv_prob_int_travel = 0.8 # ESTIMATE
		# How recent must the international travel have been to prompt further investigation?
		, int_travel_history_threshold_time = 30 # days DUMMY VALUE
		# Probability of prompting further investigation upon death
		, non_mg_inv_prob_death = 0.8 # ESTIMATE. May want to link this to age as well, i.e. only a prompt at younger ages.
		# Probability of prompting further investigation based on age of individual in ICU
		# PERHAPS USE A SLIDING SCALE - YOUNGER THE INDIVIDUAL, HIGHER THE PROBABILITY OF INVESTIGATION
		# e.g. all respiratory admissions to ICU aged 40 and under prompt investigation of pathogenesis
		# and 60+ do not, with 80% probability of investigation linear decline in probability between the two ages
		# FIND SOURCE FOR ICU ARI AGES - CAN THEN DEFINE VALUE FOR UNUSUAL AGE, e.g. X% QUANTILE OR X SDs FROM MEAN
		, non_mg_inv_prob_icu_age = vcat( fill(1,40), collect(1.0:-0.01:0.8), fill(0,40))# DUMMY VALUE TODO ADD SOURCE 
		# WOULD ALSO NEED TIMINGS AROUND TESTING AND RESULTS FOR NON-METAGENOMIC SURVEILLANCE - WILL BE LONGER - AND DIFFERENTIATE BETWEEN CURRENT ICU AND DECEASED
	);

    
    #return merge(data_dependent_params, configurable_params)
	return (default_configurable_params = configurable_params, default_P = merge(data_dependent_params, configurable_params) );
end








# Function to sample age of infectee, which is conditional on the age of the infector
function samp_infectee_age(p; contacttype, donor_age)# = :nothing, donor_age = :nothing)
	# Sample from distribution of contact ages for individual of age = x
	if isnothing( donor_age ) || (contacttype == :import)
		# i.e. this is an import and we don't know anything about the infector of this infection
		# Then we sample the infectee age from the distribution of international traveller ages
		infectee_age = Int8( wsample( INT_TRAVELLERS_AGE_SINGLE_YR[:,"age"] 
							 , INT_TRAVELLERS_AGE_SINGLE_YR[:,"intl_travel_prop_adj"] ) )
	elseif !isnothing( donor_age )
		# If there is a donor_age (infector age) for this infection then we sample the infectee age
		# according to the contact type specific contact matrix by age
		if contacttype == :F
			infectee_age = Int8( wsample( 0:(dim(p.f_contact_matrix_age) - 1) 
								  , p.f_contact_matrix_age[ donor_age + 1,:] ) )
		elseif contacttype == :G
			infectee_age = Int8( wsample( 0:(dim(p.g_contact_matrix_age) - 1) 
								  , p.g_contact_matrix_age[ donor_age + 1,:] ) )
		elseif contacttype == :H
			infectee_age = Int8( wsample( 0:(dim(P.o_contact_matrix_age) - 1)
								  , p.o_contact_matrix_age[ donor_age + 1,:] ) )
		elseif !(contacttype in CONTACTTYPES)
			error("Error in 'samp_infectee_age' function. Contact type required to return infectee age")
		end
	end
	return (infectee_age)
end
# Test
# samp_infectee_age(P; contacttype = :nothing, donor_age = 1) 
#samp_infectee_age(P; contacttype = nothing, donor_age = 1) 
#samp_infectee_age(P; contacttype = :G, donor_age = 1)
#samp_infectee_age(P; contacttype = :H, donor_age = 1)
#samp_infectee_age(P; contacttype = :F, donor_age = 1)
#samp_infectee_age(P; contacttype = :F, donor_age = nothing)

function sampdegree(p; contacttype = :nothing, age = :nothing)
	# gnegbinomr = p.gnegbinomr*p.gnegbinomp
	# gnegbinomp = p.gnegbinomp/(p.gnegbinomp + (1-p.gnegbinomp))
	# Note that element 1 of p.oorateshape (and other distribution parameter vectors) is for age = 0 years
	if age+1 > length(p.oorateshape)
		error("Error in 'sampdegree' function. No contact distribution parameters are available for individuals of this age: ", age, " years")
	end
	#Test
	#age=1
	r  = Gamma( p.oorateshape[age+1], p.ooratescale[age+1] ) |> Base.rand 
	kf = rand( NegativeBinomial(p.fnegbinomr[age+1], p.fnegbinomp[age+1]) )
	kg = rand( NegativeBinomial(p.gnegbinomr[age+1], p.gnegbinomp[age+1]) )
	if contacttype == :F # note this counts the link that transmitted infection
		kf = rand( NegativeBinomial(p.fnegbinomr[age+1]+1, p.fnegbinomp[age+1]) ) + 1
	elseif contacttype == :G
		kg = rand( NegativeBinomial(p.gnegbinomr[age+1]+1, p.gnegbinomp[age+1]) ) + 1
	elseif contacttype == :H 
		r = Gamma( p.oorateshape1[age+1], p.ooratescale1[age+1]) |> Base.rand 
	end 
	[ kf, kg, r ]	
end
#Test
#println( sampdegree(P, contacttype=:H, age = 0) )
#println( sampdegree(P, contacttype=:nothing, age = 0) ) 
#println( sampdegree(P, contacttype=:nothing, age = :nothing) )

function dayofweek(t, tinf, initialdow)
	d = Int( floor( (initialdow-1) + t-tinf ) % 7  ) + 1
	d
end
#Test
#println( dayofweek( 7, 1, 1) )

function transmissionrate(carestage, infstage, contacttype, t, tinf, tinfectious, initialdow, p, age, severity)
	dow = dayofweek(t,tinf,initialdow) 
	ρ = 1.0 # 0.250 
	# Reduced transmission while in hospital
	if carestage in (:admittedhospital,:admittedicu,:stepdown)
		ρ *= p.ρ_hosp
	end
	# Non-infectious stages
	if infstage in (:latent, :recovered, :deceased)
		ρ  *= 0.0 
	end
	# Asymptomatic individuals are less likely to transmit
	if severity.severity == :asymptomatic
		ρ *= p.ρ_asymptomatic
	end
	if contacttype == :F 
		ρ *= p.fcont
	elseif contacttype == :G
		ρ *= p.gcont
	elseif contacttype == :H
		ρ *= p.oocont
	end

	return ρ * p.dowcont[dow] * p.infectivity * pdf( Gamma(p.infectivity_shape, p.infectivity_scale), t-tinfectious )
end

function simgendist(t0, t1, p; s = 1000)
	@assert t1 >= t0 
	muts = Base.rand( 
		NegativeBinomial(s * p.μ * (t1-t0) / p.ω,  1.0 / ( 1 + p.ω ) ) # NOTE prob = 1-<paramter in ARC paper>
		)
	float( muts ) / s 
end

function simgendist(tuv, p; s = 1000)
	# Gamma(  p.μ * (tuv)/(1.0+p.ω),  1.0+p.ω ) |> Base.rand
	simgendist(0., tuv, p; s = s)
end

# Function to determine the severity of an infection
# Based on probabilities reported in Knock et al (2021), Science Translational Medicine, 13(602)
# and Saigal et al. (2025), BMC Emergency Medicine, 25:11. See 'HARISS_pathway_adjustment.jl' for calculations of some probabilities.
# Note that among patients over 65yo the probability of admission to ICU decreases with age
# but this may not be repeated in a future outbreak
# ASSUMES THAT NO-ONE DIES FROM INFECTION OUTSIDE OF HOSPITAL OR ICU (NOT IN LINE WITH COVID-19)
# Assume that ifr_by_age = p_death_icu + p_death_hosp_D + p_death_stepdown (note that this does not include deaths outside of hospital)
function sample_infectee_severity( p; age = infectee_age )
	# Begin with blank severity and fatal as false by default
	severity = nothing
	fatal = false
	# First determine if symptomatic or asymptomatic
	if rand() < (1 - p.symptomatic_prob_by_age.symptomatic_prob[ age+1 ]) # age+1 because age=0 is the first element in the vector
		severity = :asymptomatic
		fatal = false
	# If not asymptomatic then determine how severe the symptoms are
	# Severe enough for hospitalisation? 
	elseif rand() < ( p.symptomatic_ihr_by_age[age+1] )
		severity = :severe
		fatal = ( rand() < p.p_death_hosp[age+1] )
		# Severe enough for hospital and then ICU?
		# Determine if a hospitalised patient is then admitted to ICU
		if rand() < ( p.icu_by_age[age+1] )
			severity = :verysevere
			fatal = ( rand() < (p.p_death_icu[age+1] + p.p_death_stepdown[age+1])) # Probability of death if infection is very severe is sum of p of death in ICU and p of death in stepdown ward
		end

		# Split :severe infection category based on length of stay in hospital (short stay is <24h and long stay is >24h)
		if severity == :severe
			severity = StatsBase.wsample( [:severe_hosp_short_stay, :severe_hosp_long_stay], [p.prop_severe_hosp_short_stay, p.prop_severe_hosp_long_stay] )	
		end
		# Assume individuals discharged from hospital within 24h will not die from infection
		if severity == :severe_hosp_short_stay
			fatal = false
		end

	else
		# Remaining possibilities are mild, moderate_GP or moderate_ED (moderate indicates either a visit to the GP or Emergency Department but without admission to hospital or ICU)
		severity = StatsBase.wsample( [:mild, :moderate_GP, :moderate_ED], [p.prop_mild, p.prop_moderate_GP, p.prop_moderate_ED]  )
		fatal = false
	end

	infectee_severity = ( severity = severity, fatal = fatal )

	return( infectee_severity )
end
# Test
#sample_infectee_severity( P; age = 80 )

function Infection(p; pid = "0"
					, region="TLI3"
					, tinf = 0.0
					, initialdow = 1
					, contacttype = :import
					, donor::Union{Nothing,Infection} = nothing) #, age = nothing)
	# initial migration status 
	# iregion = findfirst(REGKEY.code .== region) 
	importedinfection = false
	initialregion = region
	if contacttype == :F 
		homeregion = region 
# @show homeregion
# @show keys(COMMUTEPROB[homeregion] )
# ("na" ∉ keys(COMMUTEPROB[homeregion])) && ( @bp  )
		infectee_age = samp_infectee_age( P; contacttype = contacttype
										   , donor_age = isnothing(donor) ? nothing : donor.infectee_age )
		# ONS Commuting data only for ages 16 and above. Only small numbers of commuters in the 'Aged 65y and older' category  
		iscommuter = ( infectee_age < 16 || infectee_age >= 65 ) ? false : rand() < COMMUTEPROB[homeregion]["na"] # Test # infectee_age = 10 # infectee_age = 65 # infectee_age = 64
		#iscommuter = rand() < COMMUTEPROB[homeregion]["na"] 
		prd = deepcopy( COMMUTEPROB[region] ) 
		("na" in prd.index2name) && (delete!( prd, "na" ))
		commuteregion = iscommuter ? wsample( prd.index2name, prd.data ) : region
	elseif contacttype == :import 
		importedinfection = true
		homeregion = wsample( ITL2SIZE.ITL225CD, ITL2SIZE.total_population_2022 )
		
		# Age of individual in imported infection is determined by age distribution of international travellers
		# entering the UK. 
		# Weighted sample of single year age weighted by age group and adjusted by single year in 
		# line with UK ONS population data for that age group
		infectee_age = Int8( StatsBase.wsample( INT_TRAVELLERS_AGE_SINGLE_YR[:,"age"]
				  							  , INT_TRAVELLERS_AGE_SINGLE_YR[:,"intl_travel_prop_adj"] ) )

		# ONS Commuting data only for ages 16 and above. Only small numbers of commuters in the 'Aged 65y and older' category 
		iscommuter = ( infectee_age < 16 || infectee_age >= 65 ) ? false : rand() < COMMUTEPROB[homeregion]["na"] # Test # infectee_age = 10 # infectee_age = 65 # infectee_age = 64
		#iscommuter = rand() < COMMUTEPROB[homeregion]["na"] 
		prd = deepcopy( COMMUTEPROB[region] ) 
		("na" in prd.index2name) && (delete!( prd, "na" ))
		commuteregion = iscommuter ? wsample( prd.index2name, prd.data )  : homeregion 
		contacttype = :H
	else # G or H 
		infectee_age = samp_infectee_age(P; contacttype = contacttype, donor_age = isnothing(donor) ? nothing : donor.infectee_age )
		# ONS Commuting data only for ages 16 and above. Only small numbers of commuters in the 'Aged 65y and older' category 
		iscommuter = ( infectee_age < 16 || infectee_age >= 65 ) ? false : true # Test # infectee_age = 10 # infectee_age = 65 # infectee_age = 64
		#iscommuter = true 
		prd = deepcopy( COMMUTEINPROB[region] )
		("na" in prd.index2name) && (delete!( prd, "na" ))
		homeregion  = wsample( prd.index2name, prd.data )
		commuteregion = region 
	end
	
	# transmission line list 
	H = DataFrame(pid1 = String[]
	       , pid2 = String[]
	       , timetransmission = Float64[]
	       , dayofweek = Int64[]
	       , region = String[]
	       , transmissiontype = Symbol[]
	       # , gendist0 = Float64[]
	       , degreef = Int64[]
	       , degreeg = Int64[] 
	       , oorate = Float64[] 
	       , carestage = Symbol[]
	)

	# initial state of infection 
	sampled = false 
	carestage = :undiagnosed 
	infstage = :latent # jump process will actually start with this :infectious
	isdeceased = false
	
	# initial contact network 
	flinks, glinks, hr = sampdegree(p; contacttype = contacttype , age = infectee_age )
	Kf = rand( Poisson( flinks ) )
	Kg = rand( Poisson( glinks ) )
	if contacttype == :F 
		Kf = max( Kf-1, 0)
	elseif contacttype == :G
		Kg = max( Kg-1, 0)
	end 

	# time of main events 
	tseq = Inf # time of sequencing 
	ticu = Inf # time admitted to ICU
	tgp = Inf # time visit GP 
	ted = Inf # time attend Emergency Department (ED)
	thospital = Inf # time admitted to hospital general ward
	tstepdown = Inf # time at which patient care stepped down from ICU to general ward
	tdischarge = Inf # time discharged from hospital (general ward, ICU, or stepdown after ICU)
	trecovered = Inf # time recovered - either community, discharge from hospital, ICU, or stepdown ward
	tdeceased = Inf # time of death - either in hospital, ICU, or stepdown ward

	#= Determine severity of infection, stratified by age. Probabilities from Knock et al (2021)
	Possibilities are:
	:asymptomatic (reduced transmission rate)
	:mild
	:moderate_GP (will visit GP)
	:moderate_ED (will visit Emergency Department)
	:severe_hosp_short_stay (will be admitted to hospital and stay less than 24h)
	:severe_hosp_long_stay (will be admitted to hospital and stay longer than 24h)
	:verysevere (will be admitted to ICU)
	Returns infection severity and whether infection will be fatal or not
	=# 
	severity = sample_infectee_severity( p; age = infectee_age ) #severity = StatsBase.wsample( SEVERITY, [p.propmildorasymptomatic, p.propmoderate , p.propsevere, p.propverysevere]  ) #StatsBase.wsample( SEVERITY, [P.propmild, 1-P.propmild-P.propsevere , P.propsevere]  )
	#Test
	#severity = sample_infectee_severity( P; age = 90 )
	#severity.severity
	#severity.fatal

	#cumulative transm 
	R = 0 

	# Determine length of infectious period
	gammalatent = Gamma(p.latent_shape, p.latent_scale) # median( Gamma(P.latent_shape, P.latent_scale) ) = 3.19154
	latenthazard(t) = pdf(gammalatent,t) / (1 - cdf(gammalatent,t))
	gammarecovery = Gamma(p.infectious_shape, p.infectious_scale) # median( Gamma(P.infectious_shape, P.infectious_scale) ) = 23.7
	recoveryhazard(t) = pdf(gammarecovery ,t) / (1 - cdf(gammarecovery,t)) 
	
	laglatent = rand( gammalatent )
	lagrecovery = rand( gammarecovery ) 
	tinfectious = tinf + laglatent
	tfin = laglatent + lagrecovery + tinf
	#trecovered = tfin
	if severity.fatal == false 
		trecovered = tfin
	elseif severity.fatal == true
		trecovered = Inf
		tfin = min(tdeceased,tfin)
	end 
	tspan = (tinfectious, tfin)
	
	# TODO possibly define infstage = :recovered here
	#if t >= trecovered
	#	infstage = :recovered
	#end

	# rate interval for variable rate jumps
	rint(u,p,t) = 1.0 

	# network dynamics 
	rate_gainf(u, p, t) = flinks * p.frate
	aff_gainf!(int) = begin Kf+=1 end 	
	j_gainf = ConstantRateJump(rate_gainf, aff_gainf!) 

	rate_losef(u, p, t) = Kf * p.frate 
	aff_losef!(int) = begin Kf-=1 end 
	j_losef = ConstantRateJump(rate_losef, aff_losef!) 
	
	rate_gaing(u, p, t) = glinks * p.grate
	aff_gaing!(int) = begin Kg+=1 end 
	j_gaing = ConstantRateJump(rate_gaing, aff_gaing!) 

	rate_loseg(u, p, t) = Kg * p.grate 
	aff_loseg!(int) = begin Kg-=1 end 
	j_loseg = ConstantRateJump(rate_loseg, aff_loseg!) 

	# transmissions 
	rate_transmf(u,p,t) = (region==homeregion) ? Kf*transmissionrate(carestage, infstage, :F, t, tinf, tinfectious, initialdow, p, infectee_age, severity) : 0.0
	hrate_transmf(u,p,t) = max(1.2*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p, infectee_age, severity) )
	lrate_transmf(u,p,t) = min(0.8*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p, infectee_age, severity) )
	aff_transmf!(int) = begin 
		Kf -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :F, flinks, glinks, hr, carestage)
		)
	end
	j_transmf = VariableRateJump(rate_transmf, aff_transmf!; lrate=lrate_transmf, urate=hrate_transmf, rateinterval=rint) # 
	
	rate_transmg(u,p,t) = (region==commuteregion) ? Kg*transmissionrate(carestage, infstage, :G, t,tinf,  tinfectious, initialdow, p, infectee_age, severity) : 0.0 
	hrate_transmg(u,p,t) = max(1.2*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3, tinf, tinfectious, initialdow, p, infectee_age, severity) )
	lrate_transmg(u,p,t) = min(0.8*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3,tinf,  tinfectious, initialdow, p, infectee_age, severity) )
	aff_transmg!(int) = begin 
		Kg -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :G, flinks, glinks, hr, carestage)
		)
	end
	j_transmg = VariableRateJump(rate_transmg, aff_transmg!; lrate=lrate_transmg, urate=hrate_transmg, rateinterval=rint)

	rate_transmh(u,p,t) = hr*transmissionrate(carestage, infstage, :H, t, tinf, tinfectious, initialdow, p, infectee_age, severity) 
	hrate_transmh(u,p,t) = max(1.2*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p, infectee_age, severity) )
	lrate_transmh(u,p,t) = min(0.8*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p, infectee_age, severity) )
	aff_transmh!(int) = begin 
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :H, flinks, glinks, hr, carestage)
		)
	end
	j_transmh = VariableRateJump(rate_transmh, aff_transmh!; lrate=lrate_transmh, urate=hrate_transmh, rateinterval=rint)

	## GP visit only
	rate_gp_only(u,p,t) = ((carestage==:undiagnosed) & (severity.severity in (:moderate_GP,))) ? p.gp_only_rate : 0.0  
	aff_gp_only!(int) = begin 
		carestage = :GP; 
		tgp = int.t 
	end
	j_gp_only = ConstantRateJump(rate_gp_only, aff_gp_only!)

	## GP visit before ED or hospital
	rate_gp_before_hosp(u,p,t) = ((carestage==:undiagnosed) & (severity.severity in (:moderate_ED,:severe_hosp_short_stay,:severe_hosp_long_stay,:verysevere))) ? p.gp_before_hosp_rate : 0.0  
	aff_gp_before_hosp!(int) = begin 
		carestage = :GP; 
		tgp = int.t 
	end
	j_gp_before_hosp = ConstantRateJump(rate_gp_before_hosp, aff_gp_before_hosp!)

	## Emergency Department (ED) direct visit only
	rate_ed_direct(u,p,t) = ((carestage==:undiagnosed) & (severity.severity in (:moderate_ED,))) ? p.ed_direct_rate : 0.0  
	aff_ed_direct!(int) = begin 
		#carestage = :ED; 
		ted = int.t 
		tdischarge = ted + rand( Uniform(0, p.tdischarge_ed_upper_limit) ) # Assume discharged within half a day
		carestage = :discharged # Individual is moved to this carestage because it removes the need for another jump process

	end
	j_ed_direct = ConstantRateJump(rate_ed_direct, aff_ed_direct!)

	## Emergency Department (ED) visit (but not admitted to hospital) after visiting GP
	# These patients are discharged directly from ED
	rate_ed_from_gp(u,p,t) = ((carestage==:GP) & (severity.severity in (:moderate_ED,))) ? p.ed_from_gp_rate : 0.0  
	aff_ed_from_gp!(int) = begin 
		#carestage = :ED 
		ted = int.t 
		tdischarge = ted + rand( Uniform(0, p.tdischarge_ed_upper_limit) ) # Assume discharged within half a day
		carestage = :discharged # Individual is moved to this carestage because it removes the need for another jump process
	end
	j_ed_from_gp = ConstantRateJump(rate_ed_from_gp, aff_ed_from_gp!)

	## Hospital admission direct without visit to GP
	rate_hosp_admit_direct(u,p,t) = (( (carestage in (:undiagnosed,)) & (severity.severity in (:severe_hosp_short_stay,:severe_hosp_long_stay,:verysevere)) )) ? p.hosp_admit_direct_rate : 0.0 
	aff_hosp_admit_direct!(int) = begin 
		carestage = :admittedhospital
		thospital = int.t
		if severity.severity == :severe_hosp_short_stay
			# tdischarge = thospital + random x value between 0 and 1 from the Erlang(1, 0.0935) = Gamma(1,1/0.0935)
			#tdischarge = thospital + quantile( Erlang(1, 0.0935), rand() * cdf( Erlang(1, 0.0935), 1.0 ) ) # Short stay admission to hospital is <24h (1day), which can be computed as the quantile of the Erlang distribution from Knock et al (2021) at the % using the CDF at x=1 multiplied by a random number, thereby giving a random value from the Erlang distribution between 0 and 1
			# However, the difference between this and rand(Uniform(0,1)) is not large and so the simpler form is used. The simpler form will also not need changing for different pathogens.
			tdischarge = thospital + rand( Uniform(0, p.tdischarge_hosp_short_stay_upper_limit) ) # Short stay admission to hospital is <24h (1day)
			carestage = :discharged # Individual is moved to this carestage prematurely but it stops them being available to move to :ICU or to :deceased care stages
		end
	end 
	j_hosp_admit_direct = ConstantRateJump( rate_hosp_admit_direct, aff_hosp_admit_direct! )

	## Hospital admission after visit to GP
	rate_hosp_admit_from_gp(u,p,t) = (( (carestage in (:GP,)) & (severity.severity in (:severe_hosp_short_stay,:severe_hosp_long_stay,:verysevere)) )) ? p.hosp_admit_from_gp_rate : 0.0 
	aff_hosp_admit_from_gp!(int) = begin 
		carestage = :admittedhospital  
		thospital = int.t 
		if severity.severity == :severe_hosp_short_stay
			# tdischarge = thospital + random x value between 0 and 1 from the Erlang(1, 0.0935) = Gamma(1,1/0.0935)
			#tdischarge = thospital + quantile( Erlang(1, 0.0935), rand() * cdf( Erlang(1, 0.0935), 1.0 ) ) # Short stay admission to hospital is <24h (1day), which can be computed as the quantile of the Erlang distribution from Knock et al (2021) at the % using the CDF at x=1 multiplied by a random number, thereby giving a random value from the Erlang distribution between 0 and 1
			# However, the difference between this and rand(Uniform(0,1)) is not large and so the simpler form is used. The simpler form will also not need changing for different pathogens.
			tdischarge = thospital + rand( Uniform(0, p.tdischarge_hosp_short_stay_upper_limit) ) # Short stay admission to hospital is <24h (1day)
			carestage = :discharged # Individual is moved to this carestage prematurely but it stops them being available to move to ICU or to :deceased care stage
		end
	end 
	j_hosp_admit_from_gp = ConstantRateJump( rate_hosp_admit_from_gp, aff_hosp_admit_from_gp! )

	## hospital -> discharged rate (for long-stay patients)
	rate_hosp_disch_long_stay(u,p,t) = (( (carestage in (:admittedhospital,)) & (severity.severity in (:severe_hosp_long_stay,)) & (severity.fatal == false) )) ? p.hosp_long_stay_recovery_rate : 0.0 
	aff_hosp_disch_long_stay!(int) = begin 
		# Only discharge if been in hospital for more than 24h (1day)
		if int.t > (thospital + 1)
			carestage = :discharged
			tdischarge = int.t
		end
	end 
	j_hosp_disch_long_stay = ConstantRateJump( rate_hosp_disch_long_stay, aff_hosp_disch_long_stay! )

	# This jump has been removed (possibly temporarily) because the same care stage jump and tdischarge definition is given by j_hosp_admit_from_gp and j_hosp_admit_direct
	## hospital -> discharged rate (for short-stay patients)
	#rate_hosp_disch_short_stay(u,p,t) = (( (carestage in (:admittedhospital,)) & (severity.severity in (:severe_hosp_short_stay,)) & (severity.fatal == false) )) ? p.hosp_short_stay_recovery_rate : 0.0 
	#aff_hosp_disch_short_stay!(int) = begin 
	#	carestage = :discharged
		# If time in hospital is longer than 1 day (24h) based on this jump process then force discharge time to be 1day (24h) after admission (the limit for a 'short-stay' by definition)
	#	tdischarge = int.t > (thospital + 1) ? thospital + 1 : int.t 
	#end 
	#j_hosp_disch_short_stay = ConstantRateJump( rate_hosp_disch_short_stay, aff_hosp_disch_short_stay! )

	## hospital -> death rate 
	#ratehospitaldeath(u,p,t) = (( (carestage in (:admittedhospital,)) & (severity.severity in (:severe_hosp_short_stay,:severe_hosp_long_stay)) & (severity.fatal == true) )) ? P.hosp_death_rate : 0.0 
	ratehospitaldeath(u,p,t) = (( (carestage in (:admittedhospital,)) & (severity.severity in (:severe_hosp_long_stay,)) & (severity.fatal == true) )) ? p.hosp_death_rate : 0.0 
	aff_hosp_death!(int) = begin 
		carestage = :deceased
		infstage = :deceased  
		tdeceased = int.t
		isdeceased = true
	end 
	j_hospital_death = ConstantRateJump( ratehospitaldeath, aff_hosp_death! )

	## Triage to ICU
	#rateicu(u,p,t) = (( (carestage in (:undiagnosed,:GP,:admittedhospital)) & (severity in (:verysevere,)) )) ? p.icurate : 0.0 
	#rateicu(u,p,t) = (( (carestage in (:undiagnosed,:GP,:admittedhospital)) & (severity.severity in (:verysevere,)) )) ? p.triage_icu_rate : 0.0 
	rateicu(u,p,t) = ( (carestage == :admittedhospital) & (severity.severity in (:verysevere,)) ) ? p.triage_icu_rate : 0.0 
	aff_icu!(int) = begin 
		carestage = :admittedicu  
		ticu = int.t 
		sampled = (rand() < p.p_sampled_icu ) #p.psampled )
	end 
	j_icu = ConstantRateJump( rateicu, aff_icu! )

	## ICU -> Death
	rate_icu_death(u,p,t) = (( (carestage in (:admittedicu,)) & (severity.severity in (:verysevere,)) & (severity.fatal == true ) )) ? p.icu_to_death_rate : 0.0 
	aff_icu_death!(int) = begin 
		carestage = :deceased
		infstage = :deceased
		tdeceased = int.t
		isdeceased = true 
	end 
	j_icu_death = ConstantRateJump( rate_icu_death, aff_icu_death! )

	## ICU -> Stepdown (to hospital general ward) (leading to recovery)
	rate_icu_stepdown(u,p,t) = (( (carestage in (:admittedicu,)) & (severity.severity in (:verysevere,)) )) ? p.icu_to_stepdown_leading_to_recovery_rate : 0.0 
	aff_icu_stepdown!(int) = begin 
		carestage = :stepdown  
		tstepdown = int.t 
	end 
	j_icu_stepdown = ConstantRateJump( rate_icu_stepdown, aff_icu_stepdown! )

	## Hospital general ward (ICU Stepdown) -> Discharge
	rate_stepdown_discharge(u,p,t) = (( (carestage in (:stepdown,)) & (severity.severity in (:verysevere,)) & (severity.fatal == false)) ) ? p.stepdown_to_recovery_after_icu_rate : 0.0 
	aff_stepdown_discharge!(int) = begin 
		carestage = :discharged
	end 
	j_stepdown_discharge = ConstantRateJump( rate_stepdown_discharge, aff_stepdown_discharge! )

	# migration  
	## model commuter traffic and return journeys 
	ratecommute(u,p,t) = iscommuter ? p.commuterate : 0.0 
	aff_commute!(int) = begin 
		 #@show homeregion 
		 #@show commuteregion 
		 #@show region
		 #@show setdiff( [homeregion, commuteregion], [region] ) 
		# (homeregion!=commuteregion) && (region= setdiff( [homeregion, commuteregion], [region] )[1] )
		if (homeregion!=commuteregion) & ( region in [homeregion, commuteregion] )
			region = setdiff( [homeregion, commuteregion], [region] )[1]
		end
	end 
	j_commute = ConstantRateJump( ratecommute, aff_commute! )

	rateimporthome(u,p,t) = (importedinfection & (region!=homeregion) & (region!=commuteregion)) ? p.commuterate : 0.0 
	aff_importhome!(int) = begin 
		region = homeregion 
	end 
	j_importhome = ConstantRateJump( rateimporthome, aff_importhome! )

	# simulate 
	## set to infectious; tspan starts from this point 
	infstage = :infectious
	simgenprob0 = DiscreteProblem([], tspan, p)
	## define jumps 
	jumps = [ j_gainf
		    , j_gaing
		    , j_losef
		    , j_loseg
		    , j_transmf
		    , j_transmg
		    , j_transmh 
			, j_gp_only
			, j_gp_before_hosp
			, j_ed_direct
			, j_ed_from_gp
			, j_hosp_admit_direct
			, j_hosp_admit_from_gp
			, j_hosp_disch_long_stay
			#, j_hosp_disch_short_stay # same jump currently included in j_hosp_admit_direct and j_hosp_admit_from_gp
			, j_hospital_death
		    , j_icu
			, j_icu_death
			, j_icu_stepdown
			, j_stepdown_discharge
			, j_commute
			]
	### include imported related jump (port-of-entry to home) 
	importedinfection && push!(jumps, j_importhome  )
	## jdep complete graph works, but is probably slower 
	jdep = repeat([collect(1:length(jumps))],length(jumps)) # complete graph 
	
	simgenprob1 = JumpProblem(simgenprob0, Coevolve(), jumps...; dep_graph=jdep)
	sol = solve(simgenprob1, SSAStepper())

	# evo distance   
	# if isnothing(donor) || isinf(tseq)
	# 	d = Inf 
	# elseif isinf(donor.tsequenced)
	# 	d = Inf 
	# else 
	# 	d = simgendist(abs((tseq-tinf) + (donor.tsequenced-tinf)), p)
	# end


	Infection(
				pid
				, isnothing(donor) ? missing : donor.pid 
				, H 
				, R
				, nothing #sol #nothing to save on memory 
				, tinf
				, tgp
				, ted
				, thospital
				, ticu
				, tstepdown
				, tdischarge
				, trecovered
				, tdeceased
				, contacttype
				, (flinks,glinks,hr)
				, iscommuter 
				, initialregion
				, homeregion
				, commuteregion
				, initialdow # day 1-7 when infection ocurred 
				, severity.severity
				, severity.fatal
				, isnothing(donor) ? 0 : donor.generation+1 # Determine generation 
				, isdeceased
				, isnothing(donor) ? missing : donor.infectee_age
				, isnothing(infectee_age) ? missing : infectee_age
				, importedinfection
				)
end

# Construct Infection object from row in h and donor Infection object
Infection(p, h::DataFrameRow, donor::Infection) = Infection(p; pid = h["pid2"]
							    							, region = h["region"]
															, tinf = h["timetransmission"]
															, initialdow = h["dayofweek"]
															, contacttype = h["transmissiontype"]
															, donor = donor #h["donor"] # donor
															)

# function Infection(p; pid = "0", region="TLI3", tinf = 0.0, initialdow = 1, contacttype = :nothing, donor::Union{Nothing,Infection} = nothing)

function simgeneration(p, prevgen::Array{Infection}; maxtime = Inf)
	length(prevgen) == 0 && return Array{Infection}([])
	Array{Infection}(
		[Infection(p, h, u) for u in prevgen for h in eachrow(u.H) if h["timetransmission"] < maxtime]
	)
end

function simtree(p; region="TLI3", initialtime=0.0, maxtime=30.0, maxgenerations::Int64=10, initialcontact=:H, max_cases=50000, max_cases_df, import_number, timports)
	# clustthreshold::Float64 = 0.005,
	#println( "$(max_cases) within simtree")
	#TEST
	#p=P
	#region="TLJ1"
	#initialtime=23.59195186679564
	#maxtime=90.0
	#maxgenerations=100
	#max_cases=100
	#initialcontact=:import
	
	@assert maxgenerations > 0
	@assert maxtime > initialtime

	simid = UUIDs.uuid1() |> string 

	# Define initial infection
	g = Array{Infection}([ 
	 	 					Infection(p; pid = "$(simid)-0", region = region, tinf = initialtime, contacttype = initialcontact)#, donor = nothing ) 
	# function Infection(p; pid = "0", region="TLI3", tinf = 0.0, initialdow = 1, contacttype = :nothing, donor::Union{Nothing,Infection} = nothing)
						])
	
	G = g 
	H = g[1].H 
	
	# Allocate max_cases to each import with earlier importations allowed a higher number of max cases
	# Linear - based on time from individual importation to maxtime of sim
	max_cases_per_import = floor( max_cases * (maxtime .- timports)[import_number] / sum((maxtime .- timports)))
	# Exponential - Based on estimate of cases generated by individual importation by maxtime using exponential growth
	#import_cases_est = []
	#t_import_to_maxtime = (maxtime.-timports)
	#for i in 1:length(t_import_to_maxtime)
	#	push!(import_cases_est, exponential_growth(N0 = 1.0, R = 2.0, Tg = 6.0, t = t_import_to_maxtime[i] ))
	#end
	#max_cases_per_import = floor( max_cases * (import_cases_est)[import_number] / sum(import_cases_est))

	# Add values to row of max_cases_df relating to sim for particular import
	push!(max_cases_df
		 , (simid # simid for particular importation
			, import_number # index for particular importation out of all importations
			, timports[import_number] # time of importation 
			, max_cases_df.max_cases[1] # total maximum number of cases
			, max_cases_df.n_imports[1] # total number of importations
			, max_cases_per_import #max_cases_df.max_cases_per_import[1] # maximum number of cases for particular importation, linked to the time of importation
			, "max cases NOT reached" # comment on whether the maximum was reached
			, Inf # the maximum time of infection, so can see what impact the max cases had on right censoring infections
			, Inf)) # the number of generations from the particular importation, so can see what impact the max cases had on right censoring infections
			
	for igen in 2:maxgenerations
		
		g = simgeneration(p, g, maxtime = maxtime)
		# g = [infection for infection in g if infection.tinf < maxtime]
		if length(g) > 0 
			H = vcat(H, vcat([x.H for x in g]...))
			G = vcat(G, g)
			
			# Stop simtree if max_cases reached and record comment and maximum time of infection
			if size(G,1) > max_cases_df.max_cases_per_import[size(max_cases_df,1)]
				max_cases_df[max_cases_df.simid .== simid, :comment] .= "max cases reached"
				max_cases_df[max_cases_df.simid .== simid, :max_tinf] .= maximum([(x.tinf) for x in G])
				max_cases_df[max_cases_df.simid .== simid, :generations] .= maximum([(x.generation) for x in G])
				break
			end

		else 
			break
		end
	end
	
	max_cases_df[max_cases_df.simid .== simid, :generations] .= maximum([(x.generation) for x in G])

	# removed the if statement at the end because it was stopping simforest() working after age disaggregation added
	dfargs = [(u.dpid, u.pid, u.tinf, u.contacttype, u.initialregion, u.donor_age, u.infectee_age) for u in G]# if !ismissing(u.dpid)]
	D = length(dfargs) > 0 ? 
		DataFrame(dfargs, [:donor, :recipient
							, :timetransmission, :contacttype, :region
							, :donor_age, :recipient_age]) :  
		DataFrame([:donor => nothing, :recipient => nothing
					, :timetransmission => nothing, :contacttype => nothing, :region => nothing
					, :infector_age => nothing, :infectee_age => nothing])
	
	dfargs1 = [(u.pid, u.tinf, u.tgp, u.ted, u.thospital, u.ticu, u.tstepdown, u.tdischarge, u.trecovered, u.tdeceased
				, u.severity#.severity
				, u.fatal
				, u.iscommuter, u.homeregion, u.commuteregion
				, u.generation,
				u.degree...
				, u.donor_age, u.infectee_age
				, u.importedinfection) for u in G]
	Gdf = DataFrame(dfargs1
		, [:pid, :tinf, :tgp, :ted, :thospital, :ticu, :tstepdown, :tdischarge, :trecovered, :tdeceased
			, :severity
			, :fatal
			, :iscommuter, :homeregion, :commuteregion
			, :generation
			, :F, :G, :H
			, :infector_age, :infectee_age
			, :importedinfection ]
	)
	
	H.simid .= simid
	D.simid .= simid 
	Gdf.simid .= simid 

	#isequal(mc, "$(max_cases) max cases reached") ? "Finished simtree run but $(max_cases) max cases reached" : println("Finished simtree run")
	max_cases_df[max_cases_df.simid .== simid, :max_tinf][1] < Inf ? "Finished simtree run but $(max_cases) max cases reached" : println("Finished simtree run")
	
	(
		 G = Gdf 
		, D = D 
		, infections = G 
		, H = H
		, max_cases_df = max_cases_df
	)
	
end

sampleimportregion() = begin 
	region = wsample( CAAIMPORTS.ITL225CD, CAAIMPORTS.pax_2024_per_day )
	prd = deepcopy( COMMUTEINPROB[region] )
	("na" in prd.index2name) && (delete!( prd, "na" ))
	delete!( prd, region )
	homeregion = wsample( prd.index2name, prd.data )
	(; region = homeregion , regionentry = region )
end

# simulate continuous importation 
function simforest(p; initialtime=0.0, maxtime=30.0, maxgenerations::Int64=10, initialcontact=:H, importmodel=:TDist, max_cases=50000)
#TEST
#simforest(p=NBPMscape.P; initialtime=0.0, maxtime=90.0, maxgenerations=100, initialcontact=:H, importmodel=:TDist, max_cases=1000)

	if importmodel == :Poisson # for testing  
		nimports = rand(Poisson((maxtime-initialtime)*p.importrate))
		nimports = max(1, nimports )
		timports = map(_->rand(Uniform(initialtime,maxtime)), 1:nimports)
	elseif importmodel == :TDist  # based on 1st wave covid 
		nimports = 0
		while nimports == 0 
			timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
			# use theoretical quantiles to prohibit outliers messing things up 
			tlb = quantile( TDist( p.import_t_df ), 1/p.nimports )*p.import_t_s 
			timports = timports[ (timports .> tlb) .& (timports .< (tlb+maxtime) ) ]
			nimports = length( timports )
		end
	end
	sort!(timports)
	timports .-= minimum( timports )
	importregions = map(_->sampleimportregion(),1:nimports)

	# Initialise max_cases df
	max_cases_df = DataFrame( simid = "NA"
							, import_number = Inf
							, timport = Inf
							, max_cases = max_cases
							, n_imports = nimports
							, max_cases_per_import = Inf # #floor(max_cases / nimports)
							, comment = "$(max_cases) max cases NOT reached"
							, max_tinf = Inf
							, generations = Inf)
		
	import_number = 0
	# trs = map(t->simtree(p; initialtime=t,maxtime=maxtime,maxgenerations=maxgenerations,initialcontact=:H ), timports)
	trs = map(zip(timports, importregions)) do (t, r)
		#t=timports[1]
		#r=importregions[1]
		import_number = import_number + 1
		
		simtree(p; region=r.regionentry #r.region
						, initialtime=t
						, maxtime=maxtime, maxgenerations=maxgenerations, initialcontact=:import
						, max_cases=max_cases
						, max_cases_df = max_cases_df
						, import_number
						, timports )

	end
	G = reduce( vcat, map(x-> x.G,trs))
	D = reduce( vcat, map(x-> x.D,trs))
	(; G = G, D = D, nimports = nimports, timports = timports, max_cases_df = max_cases_df[2:end,:] )
end

infectivitytoR(ν::Real; nsims = 1000) = begin
	P = merge(NBPMscape.P, (;infectivity=ν))
	infsF = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:F), 1:nsims);
	infsG = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:G), 1:nsims);
	infsH = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:H), 1:nsims);
	inftoR(inf) = begin 
		cmh = StatsBase.countmap( inf.H.transmissiontype )
		map( k-> (k in keys(cmh)) ? cmh[k] : 0, [:F, :G, :H] )
	end
	infstoR(infs) = reduce( .+, map(inftoR,infs));
	RF = infstoR(infsF) ./ nsims ;
	RG = infstoR(infsG) ./ nsims ;
	RH = infstoR(infsH) ./ nsims ;
	NGM = [ RF RG RH ]
	@show NGM
	# eigvecs(NGM)
	evs = eigvals(NGM) 
	any( map( x-> x isa Complex , evs ) ) && (return missing)
	maximum(evs)
end

sampleforest(fo::NamedTuple, p::NamedTuple) = begin 
	# icu cases 
	G = fo.G[ isfinite.(fo.G.ticu), : ]
	# sample size 
	n = rand( Binomial( size(G,1), p.p_sampled_icu )) # p.psampled ))
	# subforest 
	G1 = G[sample( 1:size(G,1), n, replace=false ), :]

	tsample = (size(G1,1)>0) ? map( g -> rand(Uniform(g.ticu,g.trecovered)) , eachrow(G1) ) : []
	G1.tsample = tsample 

	(; 
		tsample = sort( tsample  )
		, treport = sort( tsample .+ p.turnaroundtime )
		, G = G1 
		, n = n 
		, firstsample = (n>0) ? minimum(tsample) : missing
	)
end
#sampleforest(fo, psampled::Real) = begin  P = merge(NBPMscape.P,(;psampled=psampled)); sampleforest(fo,P) end
sampleforest(fo, p_sampled_icu::Real) = begin  P = merge(NBPMscape.P,(;p_sampled_icu=p_sampled_icu)); sampleforest(fo,P) end


# Function for estimating possible number of cases for individual importations so can allocate the 
# maximum case limit between multiple importations with different timport in simtree 
function exponential_growth(;N0::Float64, R::Float64, Tg::Float64, t::Float64)#::Vector{Float64})
    r = log(R) / Tg  # Growth rate derived from R and generation time
    #N_t = N0 .* exp.(r .* t)  # Element-wise exponential growth
	#N_t = N0 * exp(r * maximum(t))  # Total number of cases under exponential growth
	N_t = N0 * exp(r * t)  # Total number of cases under exponential growth
    return N_t
end

# Exp growth parameters
#N0 = 1.0 # Initial number of cases
#R = 2infectivitytoR(NBPMscape.P.infectivity; nsims = 1000) # Reproduction number
#Tg = 6.0              # Generation time in days
#t = 0:1:30            # Time vector from day 0 to day 30
#t = 0:1:maxtime            # Time vector from day 0 to day 30
#t = 0:1:200            # Time vector from day 0 to day 30
#t = (maxtime.-timports)            # Time vector from day 0 to day 30
# exponential_growth(;N0=1.0,R=2.0,Tg=6.0,t=150.0)