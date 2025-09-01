Plots.default(;lw = 4)

const MAXDURATION = 180.0

mutable struct Infection
	pid::String # infection id 
	dpid::Union{String,Missing}
	H::DataFrame # transmission linelist 
	R::Int64 # number transmissions 
	sol::Union{ODESolution,Nothing}
	tinf::Float64 # time infected 
	tgp::Float64 # time 
	thospital::Float64 # time 
	ticu::Float64 # time 
	# tsampled::Float64 # time
	# tseqdeposition::Float64 # time for sequencing + qc + bioninf 
	trecovered::Float64
	# tdeceased::Float64 # time infected individual died TODO
	contacttype::Symbol # cause of infection
	degree::Tuple{Int64,Int64,Float64} # TODO add a 4th type if change the number of contact types
	iscommuter::Bool
	initialregion::String
	homeregion::String
	commuteregion::String
	initialdayofweek::Int # day 1-7 when infection ocurred 
	severity::Symbol
	generation::Int64
	#isdeceased::Bool # TODO whether individual has died (default = FALSE). Assume all individuals that die have been admitted to ICU. 
	donor_age::Union{Int8,Missing} # Age of infector
	infectee_age::Union{Int8,Missing} # using single year age so can incorporate different age groupings for different data inputs, e.g. contact number distributions may have different age groups to the number of ICU admissions
	#importedinfection::Bool # TODO record whether individual was infected when arrived in UK - if they are admitted to ICU, this may prompt further investigation
end

#const AGEGROUPS = (:youth, :youngadult, :adult, :old, :elderly) # TODO POSSIBLY CHANGE BASED ON SC2 / NO LONGER NEEDED AS USING SINGLE YEAR AGES
const CARE  = (:undiagnosed, :GP, :admittedhospital, :admittedicu, :discharged) # :deceased) # TODO ADD DECEASED HERE??
const SEVERITY = [:mildorasymptomatic, :moderate, :severe ] # TODO POSSIBLY CHANGE BASED ON SC2
const STAGE = (:latent, :infectious, :recovered) # :deceased) # TODO ADD DECEASED HERE??

const CONTACTTYPES = (:F, :G, :H) # TODO POSSIBLY SPLIT INTO 4 (OR MORE) CONTACT TYPES: HOME, WORK/SCHOOL, TRAVEL, OTHER

# # temporary migration matrix 
# MIGMATRIX = fill( 1.0, NREGIONS, NREGIONS )
# for i in 1:NREGIONS
# 	MIGMATRIX[i,i] = 100.0 
# 	MIGMATRIX[i,:] ./= sum( MIGMATRIX[i,:] )
# end
# MIGRATES = [ sum(MIGRATES[i,:])-MIGRATES[i,i] for i in 1:NREGIONS ]

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
#=
TODO
Currently includes contact types: home, work/school, and other
Potentially change to: home, work, school, leisure, transport, other as per POLYMOD
If split work and school, may need to account for profession, i.e. teachers would have a high contact number with children 
at work everyone else would be low so would not make sense to use the average of these which would still be quite low
=#
contact_rate_dist_par_age_groups = DataFrame(
    age_group    = contact_age_groups_range
	# f is household size distribution 
	# TODO Could replace distribution derived from POLYMOD UK with that derived from ONS 
	# TODO disaggregate ONS by age
    ,fnegbinomr   = sort( filter( row -> row.contact_setting =="home", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    ,fnegbinomp   = sort( filter( row -> row.contact_setting =="home", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	# g is Work / School contact number distribution
	# POLYMOD approximate work or school contact number distribution
	# TODO disaggregate work and school because very different when disaggregated by age
    ,gnegbinomr   = sort( filter( row -> row.contact_setting =="work_school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    ,gnegbinomp   = sort( filter( row -> row.contact_setting =="work_school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	# TODO Consider splitting out School contacts in rest of model
	#,srateshape  = sort( filter( row -> row.contact_setting =="school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    #,sratescale  = sort( filter( row -> row.contact_setting =="school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	#,snegbinomr  = sort( filter( row -> row.contact_setting =="school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    #,snegbinomp  = sort( filter( row -> row.contact_setting =="school", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	# TODO Consider splitting out Work contacts in rest of model
	#,wrateshape  = sort( filter( row -> row.contact_setting =="work", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    #,wratescale  = sort( filter( row -> row.contact_setting =="work", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	#,wnegbinomr  = sort( filter( row -> row.contact_setting =="work", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    #,wnegbinomp  = sort( filter( row -> row.contact_setting =="work", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	# Other contact number distribution
	# POLYMOD approximate contact number distributions for combined settings: transport, leisure, otherplace
    ,oorateshape  = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    ,ooratescale  = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
    # excess rate distribution (also gamma) oorateshape1 = oorateshape + 1 and ooratescale1 = ooratescale
	,oorateshape1 = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape .+ 1
	,ooratescale1 = sort( filter( row -> row.contact_setting =="other_t_l_o", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	# TODO Possibly add more contact types
	# TODO Need to choose distribution
	## Transport contact number distribution (POLYMOD)
    #,trateshape  = sort( filter( row -> row.contact_setting =="transport", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    #,tratescale  = sort( filter( row -> row.contact_setting =="transport", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	#,tnegbinomr  = sort( filter( row -> row.contact_setting =="transport", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    #,tnegbinomp  = sort( filter( row -> row.contact_setting =="transport", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	## Leisure contact number distribution (POLYMOD)
    #,lrateshape  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    #,lratescale  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	#,lnegbinomr  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    #,lnegbinomp  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
	## Otherplace contact number distribution (POLYMOD)
    #,lrateshape  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_shape
    #,lratescale  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).gamma_scale
	#,lnegbinomr  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_r
    #,lnegbinomp  = sort( filter( row -> row.contact_setting =="leisure", CONTACT_DISTRIBUTIONS ), :age_group, rev = false ).nbinom_p
		
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

# TODO Note that some of these distribution parameters may not yield sensible numbers of contacts when sampled
# e.g. work/school contacts for 75plus age group
# Need way to deal with this

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

# default parameters 
global P
P = ( 
	# Relative contact rates
	# Normalised to 1 from the 10.76 contact hours in Table S2 of 
	# Danon et al (2013), Proceedings of the Royal Society B: Biological Sciences, 280(1765)
	fcont = 10.76 / 10.76   # for flinks (household) 
	, gcont = 6.71 / 10.76  # for glinks (work/school)
	, oocont = 8.09 / 10.76 # for other
	# TODO POSSIBLY ADD 4TH CONTACT TYPE (TRAVEL) AS PER WARWICK UNI SOCIAL CONTACT SURVEY (Danon et al (2013))
	#, tcont = 0.45 / 10.76  # for travel links

	# Scale adjustment for contact rates by day of the week (Sun-Sat)
	# Sourced from POLYMOD - Mossong et al (2008), PLOS Medicine, 5(3)
	# TODO COULD DISAGGREGATE THIS BY CONTACT TYPE, E.G. SCHOOL/WORK WOULD BE VERY DIFFERENT FROM LEISURE
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
	

	# TODO introduce population variation in infectiousness 
	, infectivity = 1.25 # scales transmission rate 
	, infectivity_shape = 2.2 * 0.75 # TODO 
	, infectivity_scale = 2.5 * 0.75
	# , infectivity_shape = 1.72 # Jones & Drosten, Science 2021 
	# , infectivity_scale = 5.95

	# TODO STRATIFY BY AGE?
	# , latent_shape = 4.24 # Galmiche Lancet Microbe 2023, mean ~ 5d
	# , latent_scale = 1.08 
	, latent_shape = 3.26 # Zhao 2021 
	, latent_scale = 0.979
	, infectious_shape = 8.16 # Verity 2020 , mean 24d 
	, infectious_scale = 3.03 

	, ρ = 0.250 #  transmission reduction 
	
	, frate = 0.0 # rate of gaining & losing flinks
	, grate = 1/30.0 # rate of gaining and losing
	# TODO DO WE NEED TO ADD RATES FOR OTHER AND TRAVEL?
	#, orate = TODO
	#, trate = TODO

	# Contact rate distributions stratified by age of individual and contact type
	# Based on POLYMOD UK - Mossong et al (2008), PLOS Medicine, 5(3)
	# Note that element 1 is for age = 0 years
	, fnegbinomr   = contact_rate_dist_par_all_ages[:fnegbinomr]
	, fnegbinomp   = contact_rate_dist_par_all_ages[:fnegbinomp]
	, gnegbinomr   = contact_rate_dist_par_all_ages[:gnegbinomr]
	, gnegbinomp   = contact_rate_dist_par_all_ages[:gnegbinomp]
	#, trateshape  = contact_rate_dist_par_all_ages[:trateshape]
	#, tratescale  = contact_rate_dist_par_all_ages[:tratescale]
	, oorateshape  = contact_rate_dist_par_all_ages[:oorateshape]
	, ooratescale  = contact_rate_dist_par_all_ages[:ooratescale]
	, oorateshape1 = contact_rate_dist_par_all_ages[:oorateshape1]
	, ooratescale1 = contact_rate_dist_par_all_ages[:ooratescale1]

	# Contact rate matrix by age
	# Based on POLYMOD UK - Mossong et al (2008), PLOS Medicine, 5(3)
	# TODO COULD SPLIT INTO MORE CONTACT TYPES
	, f_contact_matrix_age = CONTACT_MATRIX_HOME_SINGLE_YEAR
	, g_contact_matrix_age = CONTACT_MATRIX_SCHOOL_WORK_SINGLE_YEAR
	, o_contact_matrix_age = CONTACT_MATRIX_OTHER_SINGLE_YEAR

	# Estimated uniform delay from sampling to sequencing+bioinformatics+database
	# TODO SEARCH FOR SOURCE
	, lagseqdblb  = 3 # lower bound TODO THIS DOESN'T SEEM TO BE USED ANYWHERE
	, lagsseqdbub = 7 # upper bound TODO THIS DOESN'T SEEM TO BE USED ANYWHERE
	
	# TODO STRATIFY BY AGE GROUP - SEE KNOCK ET AL (2021)
	, propmild = 0.60 
	, propsevere = 0.05

	# TODO STRATIFY BY AGE GROUP
	, gprate = 1/3 
	, hospadmitrate = 1/4 # Docherty 2020 
	, icurate = 1/2.5 # Knock 2021
	#, hosp_disch_rate # TODO - SEE KNOCK ET AL (2021)
	#, icu_disch_rate # TODO - SEE KNOCK ET AL (2021)
	#, death_rate =  TODO ADD RATE FOR DEATHS - POSSIBLY DIFFERENT RATES DEPENDING ON CARE STAGE - SEE KNOCK ET AL (2021)

	, psampled = .05  # prop sampled from icu 
	, turnaroundtime = 3 # days TODO HOW DOES THIS LINK TO lagseqdblb and lagsseqdbub?

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

	## Assumptions regarding current non-metagenomic surveillance
	# Proportion of ICU admissions with history of international travel that would be 
	# further investigated (with non-metagenomic surveillance)
	, non_mg_inv_prob_int_travel = 0.8 # ESTIMATE TODO ADD SOURCE. HOW RECENT MUST INT'L TRAVEL BE. NEED DATA ON INT'L TRAVEL FOR POPULATION.
	# How recent must the international travel have been to prompt further investigation?
	, int_travel_history_threshold_time = 30 # days DUMMY VALUE TODO ADD SOURCE
	# Probability of prompting further investigation upon death
	# TODO ANY DEATH OR ONLY IF COMBINED WITH YOUNGER AGE AND/OR INT'L TRAVEL
	, non_mg_inv_prob_death = 0.8 # ESTIMATE TODO ADD SOURCE
	# Probability of prompting further investigation based on age of individual in ICU
	# TODO PERHAPS USE A SLIDING SCALE - YOUNGER THE INDIVIDUAL, HIGHER THE PROBABILITY OF INVESTIGATION
	# E.G. all respiratory admissions to ICU aged 40 and under prompt investigation of pathogenesis
	# and 60+ do not, with 80% probability of investigation linear decline in probability between the two ages
	# FIND SOURCE FOR ICU ARI AGES - CAN THEN DEFINE VALUE FOR UNUSUAL AGE, E.G. X% QUANTILE OR X SDs FROM MEAN
	, non_mg_inv_prob_icu_age = vcat( fill(1,40), collect(1.0:-0.01:0.8), fill(0,40))# DUMMY VALUE TODO ADD SOURCE 
	# TODO NEED TIMINGS AROUND TESTING AND RESULTS FOR NON-METAGENOMIC SURVEILLANCE - WILL BE LONGER - AND DIFFERENTIATE BETWEEN CURRENT ICU AND DECEASED

)

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
		elseif isnothing( contacttype )
			infectee_age = nothing #Int8( usample )
			#TODO DO WE WANT TO RETURN SOMETHING IF CONTACT TYPE NOT DEFINED
			# E.G. SAMPLE FROM CONTACT MATRIX INCLUDING ALL CONTACT TYPES? 
		end
	end
	return (infectee_age)
end
# Test
# samp_infectee_age(P; contacttype = :nothing, donor_age = 1) # TODO CURRENTLY RETURNS AN ERROR
#samp_infectee_age(P; contacttype = :G, donor_age = 1)
#samp_infectee_age(P; contacttype = :H, donor_age = 1)
#samp_infectee_age(P; contacttype = :F, donor_age = 1)
#samp_infectee_age(P; contacttype = :F, donor_age = nothing)

function sampdegree(p; contacttype = :nothing, age = :nothing)
	# gnegbinomr = p.gnegbinomr*p.gnegbinomp
	# gnegbinomp = p.gnegbinomp/(p.gnegbinomp + (1-p.gnegbinomp))
	# Note that element 1 of p.oorateshape (and other distribution parameter vectors) is for age = 0 years
	if age+1 > length(p.oorateshape)
		println("No contact distribution parameters are available for individuals of this age: ", age, " years")
		return
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

function dayofweek(t, tinf, initialdow)
	d = Int( floor( (initialdow-1) + t-tinf ) % 7  ) + 1
	d
end
#Test
#println( dayofweek( 7, 1, 1) )

# TODO POSSIBLY STRATIFY BY AGE GROUP / SINGLE YEAR AGE
function transmissionrate(carestage, infstage, contacttype, t, tinf, tinfectious, initialdow, p, age)
	dow = dayofweek(t,tinf,initialdow) 
	ρ = 1.0 # 0.250 
	if carestage in (:admittedhospital,:admittedicu)
		ρ *= p.ρ
	end
	if infstage in (:latent, :recovered)
		ρ  *= 0.0 
	end
	if contacttype == :F 
		ρ *= p.fcont
	elseif contacttype == :G
		ρ *= p.gcont
	elseif contacttype == :H
		ρ *= p.oocont
	end

	return ρ * p.dowcont[dow] * p.infectivity * pdf( Gamma(p.infectivity_shape, p.infectivity_scale), t-tinfectious )
	# ρ * p.dowcont[dow] * p.infectivity * pdf( Gamma(p.infectivity_shape, p.infectivity_scale), t-tinf) # TODO clarify if this is docked to tinf or tinfectious
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
		# ONS Commuting data only for ages 16 and above 
		# TODO COULD TRY TO FURTHER STRATIFY COMMUTING BY AGE IF REQUIRED
		iscommuter = infectee_age < 16 ? false : rand() < COMMUTEPROB[homeregion]["na"] 
		#iscommuter = rand() < COMMUTEPROB[homeregion]["na"] 
		prd = deepcopy( COMMUTEPROB[region] ) 
		("na" in prd.index2name) && (delete!( prd, "na" ))
		commuteregion = iscommuter ? wsample( prd.index2name, prd.data ) : region
	elseif contacttype == :import 
		importedinfection = true
		homeregion = wsample( ITL2SIZE.ITL225CD, ITL2SIZE.total_population_2022 ) # TODO COULD STRATIFY BY AGE BUT PERHAPS NOT WORTHWHILE
		
		# Age of individual in imported infection is determined by age distribution of international travellers
		# entering the UK. Two possible methods: TODO CHOOSE ONE
		# (1) Weighted sample of single year age weighted by age group and adjusted by single year in 
		# line with UK ONS population data for that age group
		infectee_age = Int8( StatsBase.wsample( INT_TRAVELLERS_AGE_SINGLE_YR[:,"age"]
				  							  , INT_TRAVELLERS_AGE_SINGLE_YR[:,"intl_travel_prop_adj"] ) )
		# OR
		# (2) Weighted sample of age group... 
		#age_group = StatsBase.wsample( INT_TRAVELLERS_AGE_GROUP[:,"Age"], INT_TRAVELLERS_AGE_GROUP[:,"percent"] ) 
		# ... then uniform sample of single year age from age group range
		#if age_group == "65 & over"
		#	age = Int8( rand(65:100) )
		#else
		#	age = Int8( rand( parse_range( age_group ) ) )# Returns error for oldest age group "65 & over"
		#end

		# ONS Commuting data only for ages 16 and above 
		# TODO COULD TRY TO FURTHER STRATIFY COMMUTING BY AGE IF REQUIRED (categories in ONS data: 16-24, 25-34, 35-44, 45-54, 55-64, 65+)
		iscommuter = infectee_age < 16 ? false : rand() < COMMUTEPROB[homeregion]["na"] 
		#iscommuter = rand() < COMMUTEPROB[homeregion]["na"] 
		prd = deepcopy( COMMUTEPROB[region] ) 
		("na" in prd.index2name) && (delete!( prd, "na" ))
		commuteregion = iscommuter ? wsample( prd.index2name, prd.data )  : homeregion 
		contacttype = :H
	else # G or H 
		infectee_age = samp_infectee_age(P; contacttype = contacttype, donor_age = isnothing(donor) ? nothing : donor.infectee_age )
		# ONS Commuting data only for ages 16 and above 
		# TODO COULD TRY TO FURTHER STRATIFY COMMUTING BY AGE IF REQUIRED
		iscommuter = infectee_age < 16 ? false : true 
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
		   #, donor_age = Int8[]    # Caused issues when added here, so only recorded in D and G
		   #, infectee_age = Int8[] # Caused issues when added here, so only recorded in D and G
	)

	# initial state of infection 
	sampled = false 
	carestage = :undiagnosed 
	infstage = :latent # jump process will actually start with this :infectious
	isdeceased = false # TODO NEED TO ADD DECEASED COMPARTMENT IN MODEL
	#isrecovered = false # TODO IS THIS REQUIRED?
	
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
	ticu = Inf # icu admit 
	tgp = Inf # gpu attend 
	thospital = Inf # hospital admit 
	tsampled = Inf # sampled 
	tseqdeposition = Inf # sequence db 
	trecovered = Inf # recovered/deceased/removed # TODO ONCE DECEASED COMPARTMENT ADDED, REMOVE DECEASED
	#tdeceased = Inf # deceased TODO ADD ONCE DECEASED COMPARTMENT INCLUDED

	# TODO AGE STRATIFY SEE KNOCK ET AL (2021)
	severity = StatsBase.wsample( SEVERITY, [p.propmild, 1-p.propmild-p.propsevere , p.propsevere]  )
	
	#cumulative transm 
	R = 0 

	gammalatent = Gamma(p.latent_shape, p.latent_scale)
	latenthazard(t) = pdf(gammalatent,t) / (1 - cdf(gammalatent,t))
	gammarecovery = Gamma(p.infectious_shape, p.infectious_scale)
	recoveryhazard(t) = pdf(gammarecovery ,t) / (1 - cdf(gammarecovery,t)) 
	# TODO ADD ANOTHER FUNCTION HERE ONCE DECEASED COMPARTMENT ADDED?

	laglatent = rand( gammalatent )
	lagrecovery = rand( gammarecovery ) 
	tinfectious = tinf + laglatent
	tfin = laglatent + lagrecovery + tinf
	trecovered = tfin 
	tspan = (tinfectious, tfin)
	
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
	rate_transmf(u,p,t) = (region==homeregion) ? Kf*transmissionrate(carestage, infstage, :F, t,tinf,  tinfectious, initialdow, p, infectee_age) : 0.0 
	hrate_transmf(u,p,t) = max(1.2*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p, infectee_age) )
	lrate_transmf(u,p,t) = min(0.8*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p, infectee_age) )
	aff_transmf!(int) = begin 
		Kf -= 1
		R += 1
		nextpid = pid * ".$(R)"
		# TODO this method is not correctly recording the donor/infector and recipient/infectee ages - it is one step behind
		#donor_age = isnothing(donor) ? -1 : donor.infectee_age
		#infectee_age = isnothing(infectee_age) ? -2 : infectee_age
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :F, flinks, glinks, hr, carestage)#, donor_age, infectee_age)
		)
	end
	j_transmf = VariableRateJump(rate_transmf, aff_transmf!; lrate=lrate_transmf, urate=hrate_transmf, rateinterval=rint) # 
	
	rate_transmg(u,p,t) = (region==commuteregion) ? Kg*transmissionrate(carestage, infstage, :G, t,tinf,  tinfectious, initialdow, p, infectee_age) : 0.0 
	hrate_transmg(u,p,t) = max(1.2*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3, tinf, tinfectious, initialdow, p, infectee_age) )
	lrate_transmg(u,p,t) = min(0.8*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3,tinf,  tinfectious, initialdow, p, infectee_age) )
	aff_transmg!(int) = begin 
		Kg -= 1
		R += 1
		nextpid = pid * ".$(R)"
		# TODO this method is not correctly recording the donor/infector and recipient/infectee ages - it is one step behind
		#donor_age = isnothing(donor) ? -1 : donor.infectee_age
		#infectee_age = isnothing(infectee_age) ? -2 : infectee_age
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :G, flinks, glinks, hr, carestage)#,  donor_age, infectee_age)
		)
	end
	j_transmg = VariableRateJump(rate_transmg, aff_transmg!; lrate=lrate_transmg, urate=hrate_transmg, rateinterval=rint)

	rate_transmh(u,p,t) = hr*transmissionrate(carestage, infstage, :H, t, tinf, tinfectious, initialdow, p, infectee_age) 
	hrate_transmh(u,p,t) = max(1.2*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p, infectee_age) )
	lrate_transmh(u,p,t) = min(0.8*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p, infectee_age) )
	aff_transmh!(int) = begin 
		R += 1
		nextpid = pid * ".$(R)"
		# TODO this method is not correctly recording the donor/infector and recipient/infectee ages - it is one step behind
		#donor_age = isnothing(donor) ? -1 : donor.infectee_age
		#infectee_age = isnothing(infectee_age) ? -2 : infectee_age
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :H, flinks, glinks, hr, carestage)#,  donor_age, infectee_age)
		)
	end
	j_transmh = VariableRateJump(rate_transmh, aff_transmh!; lrate=lrate_transmh, urate=hrate_transmh, rateinterval=rint)


	# care pathway  
	# const CARE  = (:undiagnosed, :GP, :admittedhospital, :admittedicu, :discharged) 
	# const SEVERITY = (:mildorasymptomatic, :moderate, :severe )
	
	# tseq = Inf # time of sequencing 
	# ticu = Inf # icu admit 
	# tgp = Inf # gpu attend 
	# thospital = Inf # hospital admit 
	# tsampled = Inf # sampled 
	# tseqdeposition = Inf # sequence db 
	# trecovered = Inf # recovered/deceased/removed

	## gp 
	rategp(u,p,t) = ((carestage==:undiagnosed) & (severity in (:moderate,:severe))) ? p.gprate : 0.0 
	aff_gp!(int) = begin 
		carestage = :GP; 
		tgp = int.t 
	end
	j_gp = ConstantRateJump(rategp, aff_gp!)

	## hospital
	ratehospital(u,p,t) = (( (carestage in (:undiagnosed,:GP)) & (severity in (:severe,)) )) ? p.hospadmitrate : 0.0 
	aff_hosp!(int) = begin 
		carestage = :admittedhospital  
		thospital = int.t 
	end 
	j_hospital = ConstantRateJump( ratehospital, aff_hosp! )

	##TODO add hospital -> discharged rate 
	#ratehospitaldischarge(u,p,t) = (( (carestage in (:admittedhospital)) & (severity in (:severe,)) )) ? p.hosp_disch_rate : 0.0 
	#aff_hosp_disch!(int) = begin 
	#	carestage = :dischargedhospital  
	#	thospital = int.t # TODO
	#end 
	#j_hospital_discharge = ConstantRateJump( ratehospitaldischarge, aff_hosp_disch! )


	## icu
	rateicu(u,p,t) = (( (carestage in (:undiagnosed,:GP,:admittedhospital)) & (severity in (:severe,)) )) ? p.icurate : 0.0 
	aff_icu!(int) = begin 
		carestage = :icu  
		ticu = int.t 
		sampled = (rand() < p.psampled )
	end 
	j_icu = ConstantRateJump( rateicu, aff_icu! )

	## TODO icu discharge to hospital ward 
	#rate_icu_discharge(u,p,t) = (( (carestage in (:icu)) & (severity in (:severe,)) )) ? p.icu_disch_rate : 0.0 
	#aff_icu_disch!(int) = begin 
	#	carestage = :hospital_stepdown  
	#	ticu = int.t 
	#end 
	#j_icu_discharge = ConstantRateJump( rate_icu_discharge, aff_icu_disch! )

	# # sampling 
	# ratesample(u,p,t) = (carestage == :icu & sampled=false) ? p.samplerate : 0.0 
	# aff_sample!(int) = begin 
	# 	sampled = true 
	# 	tsampled = int.t 
	# end
	# j_sample = ConstantRateJump( ratesample, aff_sample! )
	

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

	rateimporthome(u,p,t) = (importedinfection & (region!=homeregion)) ? p.commuterate : 0.0 
	aff_importhome!(int) = begin 
		region = homeregion 
	end 
	j_importhome = ConstantRateJump( rateimporthome, aff_importhome! )

	## TODO ADD SET OF RATES AND AFFECTS FOR DEATHS
	# j_death = TODO

	## TODO long distance occasional travel 
	
	# ratemigration(u,p,t) = MIGRATES[region] 
	# aff_migration!(int) = begin 
	# 	w = MIGMATRIX[region, :] 
	# 	w[region] = 0.0 
	# 	region = StatsBase.wsample( 1:NREGIONS, w )
	# end
	# j_migration = ConstantRateJump( ratemigration, aff_migration! )

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
		    , j_gp 
		    , j_hospital
			#, j_hospital_discharge TODO
		    , j_icu
			#, j_icu_discharge TODO THIS IS REQUIRED BECAUSE SOME INDIVIDUALS STEPDOWN FROM ICU TO GENERAL WARD BEFORE DYING AND WE WANT TO USE DEATH AS A PROMPT FOR SAMPLING
		    #, j_migration
		    , j_commute
		   #, j_death TODO
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
				, thospital
				, ticu 
				# , tsampled 
				# , tseqdeposition 
				, trecovered
				#, tdeceased # TODO
				, contacttype
				, (flinks,glinks,hr)
				, iscommuter 
				, initialregion
				, homeregion
				, commuteregion
				, initialdow # day 1-7 when infection ocurred 
				, severity
				, isnothing(donor) ? 0 : donor.generation+1 # Determine generation 
				#, isdeceased # TODO
				, isnothing(donor) ? missing : donor.infectee_age
				, isnothing(infectee_age) ? missing : infectee_age
				#, importedinfection # TODO
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

# TODO IN SIMTREE NEED TO RECORD IF INFECTION IS AN IMPORTATION OR NOT BECAUSE IMPORTATION MAY PROMPT
# INVESTIGATION UNDER NON-METAGENOMIC SURVEILLANCE X% OF THE TIME
function simtree(p; region="TLI3", initialtime=0.0, maxtime=30.0, maxgenerations::Int64=10, initialcontact=:H)
	# clustthreshold::Float64 = 0.005,
	
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
	
	for igen in 2:maxgenerations
		g = simgeneration(p, g, maxtime = maxtime)
		# g = [infection for infection in g if infection.tinf < maxtime]
		if length(g) > 0 
			H = vcat(H, vcat([x.H for x in g]...))
			G = vcat(G, g)
		else 
			break
		end
	end
	
	# removed the if statement at the end because it was stopping simforest() working after age disaggregation added
	dfargs = [(u.dpid, u.pid, u.tinf, u.contacttype, u.initialregion, u.donor_age, u.infectee_age) for u in G]# if !ismissing(u.dpid)]
	D = length(dfargs) > 0 ? 
		DataFrame(dfargs, [:donor, :recipient
							, :timetransmission, :contacttype, :region
							, :donor_age, :recipient_age]) :  
		DataFrame([:donor => nothing, :recipient => nothing
					, :timetransmission => nothing, :contacttype => nothing, :region => nothing
					, :infector_age => nothing, :infectee_age => nothing])
	
	dfargs1 = [(u.pid, u.tinf, u.tgp, u.thospital, u.ticu, u.trecovered
				, u.severity
				, u.iscommuter, u.homeregion, u.commuteregion
				, u.generation,
				u.degree...
				, u.donor_age, u.infectee_age) for u in G]
	Gdf = DataFrame(dfargs1
		, [:pid, :tinf, :tgp, :thospital, :ticu, :trecovered
			, :severity
			, :iscommuter, :homeregion, :commuteregion
			, :generation
			, :F, :G, :H
			, :infector_age, :infectee_age ]
	)
	
	H.simid .= simid 
	D.simid .= simid 
	Gdf.simid .= simid 

	(
		 G = Gdf 
		, D = D 
		, infections = G 
		, H = H 
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
function simforest(p; initialtime=0.0, maxtime=30.0, maxgenerations::Int64=10, initialcontact=:H, importmodel=:TDist)
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
	# trs = map(t->simtree(p; initialtime=t,maxtime=maxtime,maxgenerations=maxgenerations,initialcontact=:H ), timports)
	trs = map(zip(timports, importregions)) do (t, r)
		simtree(p; region=r.regionentry #r.region
					, initialtime=t
					, maxtime=maxtime, maxgenerations=maxgenerations, initialcontact=:import)
	end
	G = reduce( vcat, map(x-> x.G,trs))
	D = reduce( vcat, map(x-> x.D,trs))
	(; G = G, D = D, nimports = nimports, timports = timports  )
end

infectivitytoR(ν::Real; nsims = 1000) = begin
	P = merge(NBPMscape.P, (;infectivity=ν))
	infsF = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:F), 1:nsims); # TODO add donor = nothing?
	infsG = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:G), 1:nsims); # TODO add donor = nothing?
	infsH = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:H), 1:nsims); # TODO add donor = nothing?
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
	n = rand( Binomial( size(G,1), p.psampled ))
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
sampleforest(fo, psampled::Real) = begin  P = merge(NBPMscape.P,(;psampled=psampled)); sampleforest(fo,P) end 
