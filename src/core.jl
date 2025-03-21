default(;lw = 4)

const MAXDURATION = 180.0

mutable struct Infection
	pid::String # infection id 
	dpid::Union{String,Missing}
	H::DataFrame # transmission linelist 
	R::Int64 # number transmissions 
	sol::ODESolution
	tinf::Float64 # time infected 
	tgp::Float64 # time 
	thospital::Float64 # time 
	ticu::Float64 # time 
	# tsampled::Float64 # time
	# tseqdeposition::Float64 # time for sequencing + qc + bioninf 
	trecovered::Float64 
	contacttype::Symbol # cause of infection
	degree::Tuple{Int64,Int64,Float64}
	iscommuter::Bool
	initialregion::String
	homeregion::String
	commuteregion::String
	initialdayofweek::Int # day 1-7 when infection ocurred 
	generation::Int64
end

const AGEGROUPS = (:youth, :youngadult, :adult, :old, :elderly) 
const CARE  = (:undiagnosed, :GP, :admittedhospital, :admittedicu, :discharged) 
const SEVERITY = [:mildorasymptomatic, :moderate, :severe ]
const STAGE = (:latent, :infectious, :recovered)

const CONTACTTYPES = (:F, :G, :H)

const NREGIONS = 50
const REGIONS = ( :TODO )

# # temporary migration matrix 
# MIGMATRIX = fill( 1.0, NREGIONS, NREGIONS )
# for i in 1:NREGIONS
# 	MIGMATRIX[i,i] = 100.0 
# 	MIGMATRIX[i,:] ./= sum( MIGMATRIX[i,:] )
# end
# MIGRATES = [ sum(MIGRATES[i,:])-MIGRATES[i,i] for i in 1:NREGIONS ]


# default parameters 
P = ( 
	fcont = 1.0 # relative contact rate for flinks  (household) # TODO get from contact tracing studies 
	, gcont = .50 #.25 # contact rate for glinks (work)
	, oocont = .50 #.05 # contact rate stranger 

	, dowcont =(0.1043502
		,0.1402675
		,0.1735913
		,0.1437642
		,0.1596205
		,0.1445298
		,0.1338766
	) #(.2, .75, 1.0, 1.0, 1.0, .8, .3) # day of week scale Sun-Sat  (polymod )
	, infectivity = 6.0 # scales transmission rate 
	, infectivity_shape = 2.2 # Metcalf Nat Comm 2021 
	, infectivity_scale = 2.5 

	, latent_shape = 2.3 # Casey 2021 ?
	, latent_scale = 1.2
	, infectious_shape = 2.1 # Verity 2020 , mean 22 d ?
	, infectious_scale = 10.57 

	, ρ = 0.250 #  transmission reduction 
	
	, frate = 0.0 # rate of gaining & losing flinks
	, grate = 1/30.0 # rate of gaining and losing

	, fnegbinomr = 2.38e6# ONS # TODO should replace with Poisson #  8.1 Approximate household size distribution  
	, fnegbinomp =  0.999999 # 0.86 
	, gnegbinomr =  1.44 # polymod # 3 # Approximate workplace size distribution 
	, gnegbinomp = 0.1366 # 0.25 
	, oorateshape = 1.42 # polymod  # Approximate other contacts (e.g. public transport) 
	, ooratescale = 6.27746 # 
	, oorateshape1 = 2.42 # excess rate distribution (also gamma)
	, ooratescale1 = 6.27746
	, lagseqdblb = 3 # Uniform delay from sampling to sequencing+bioinformatics+database
	, lagsseqdbub = 7
	
	, propmild = 0.60 
	, propsevere = 0.05

	, gprate = 1/3 
	, hospadmitrate = 1/3
	, icurate = 1/5 
	, psampled = .05  # prop sampled form icu 

	, commuterate = 2.0

	, importrate = .5 # if using constant rate 
	, nimports = 1000 
	, import_t_df = 2.48 # personal analysis of lineages studied in Volz et al. Cell 2020 
	, import_t_s = 8.75
 #      df           m           s    
 #   2.482146   -1.074888    8.749736 
 # ( 1.226155) ( 1.608095) ( 1.978287)

	, μ = 0.001 # mean clock rate -- additive relaxed clock
	, ω = 0.5 # variance inflation
)


function sampdegree(p; contacttype = :nothing)
	# gnegbinomr = p.gnegbinomr*p.gnegbinomp
	# gnegbinomp = p.gnegbinomp/(p.gnegbinomp + (1-p.gnegbinomp))
	r = Gamma( p.oorateshape, p.ooratescale ) |> Base.rand 
	kf = rand( NegativeBinomial(p.fnegbinomr, p.fnegbinomp) )
	kg = rand( NegativeBinomial(p.gnegbinomr, p.gnegbinomp) )
	if contacttype == :F # note this counts the link that transmitted infection
		kf = rand( NegativeBinomial(p.fnegbinomr+1, p.fnegbinomp) ) + 1 
	elseif contacttype == :G
		# kg = rand( NegativeBinomial(gnegbinomr+1, gnegbinomp) ) + 1 
		kg = rand( NegativeBinomial(p.gnegbinomr+1, p.gnegbinomp) ) + 1 
	elseif contacttype == :H 
		r = Gamma( p.oorateshape1, p.ooratescale1) |> Base.rand 
	end 
	[ kf, kg, r ]	
end

function dayofweek(t, tinf, initialdow)
	d = Int( floor( (initialdow-1) + t-tinf ) % 7  ) + 1
	d
end

function transmissionrate(carestage, infstage, contacttype, t, tinf, tinfectious, initialdow, p) 
	dow = dayofweek(t,tinf,initialdow) 
	ρ = 0.250 
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

	ρ * p.dowcont[dow] * p.infectivity * pdf( Gamma(p.infectivity_shape, p.infectivity_scale), t-tinfectious )
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

function Infection(p; pid = "0", region="TLI3", tinf = 0.0, initialdow = 1, contacttype = :nothing, donor::Union{Nothing,Infection} = nothing)
	# initial migration status 
	# iregion = findfirst(REGKEY.code .== region) 
	initialregion = region
	if contacttype == :F 
		homeregion = region 
# @show homeregion
# @show keys(COMMUTEPROB[homeregion] )
# ("na" ∉ keys(COMMUTEPROB[homeregion])) && ( @bp  )
		iscommuter = rand() < COMMUTEPROB[homeregion]["na"] 
		prd = deepcopy( COMMUTEPROB[region] ) 
		("na" in prd.index2name) && (delete!( prd, "na" ))
		commuteregion = iscommuter ? wsample( prd.index2name, prd.data )  : region # TODO 
	else 
		iscommuter = true 
		prd = deepcopy( COMMUTEINPROB[region] )
		("na" in prd.index2name) && (delete!( prd, "na" ))
		homeregion  = wsample( prd.index2name, prd.data )
		commuteregion = region 
	end
	# TODO elif contact is H vs G 
	
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

	# initial contact network 
	flinks, glinks, hr = sampdegree(p; contacttype = contacttype)
	Kf = rand( Poisson(flinks) )
	Kg = rand( Poisson( glinks ))
	if contacttype == :F 
		Kf = max( Kf-1,0)
	elseif contacttype == :G
		Kg = max(Kg-1, 0)
	end 

	# time of main events 
	tseq = Inf # time of sequencing 
	ticu = Inf # icu admit 
	tgp = Inf # gpu attend 
	thospital = Inf # hospital admit 
	tsampled = Inf # sampled 
	tseqdeposition = Inf # sequence db 
	trecovered = Inf # recovered/deceased/removed

	# initial state of infection 
	sampled = false 
	carestage = :undiagnosed 
	infstage = :latent # jump process will actually start with this :infectious
	agegroup = :adult # TODO sample 
	severity = StatsBase.wsample( SEVERITY, [p.propmild, 1-p.propmild-p.propsevere , p.propsevere]  )

	#cumulative transm 
	R = 0 

	gammalatent = Gamma(p.latent_shape, p.latent_scale)
	latenthazard(t) = pdf(gammalatent,t) / (1 - cdf(gammalatent,t))
	gammarecovery = Gamma(p.infectious_shape, p.infectious_scale)
	recoveryhazard(t) = pdf(gammarecovery ,t) / (1 - cdf(gammarecovery,t)) 

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
	rate_transmf(u,p,t) = (region==homeregion) ? Kf*transmissionrate(carestage, infstage, :F, t,tinf,  tinfectious, initialdow, p) : 0.0 
	hrate_transmf(u,p,t) = max(1.2*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p) )
	lrate_transmf(u,p,t) = min(0.8*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3,tinf,  tinfectious, initialdow, p) )
	aff_transmf!(int) = begin 
		Kf -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :F, flinks, glinks, hr, carestage)
		)
	end
	j_transmf = VariableRateJump(rate_transmf, aff_transmf!; lrate=lrate_transmf, urate=hrate_transmf, rateinterval=rint) # 
	
	rate_transmg(u,p,t) = (region==commuteregion) ? Kg*transmissionrate(carestage, infstage, :G, t,tinf,  tinfectious, initialdow, p) : 0.0 
	hrate_transmg(u,p,t) = max(1.2*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3, tinf, tinfectious, initialdow, p) )
	lrate_transmg(u,p,t) = min(0.8*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3,tinf,  tinfectious, initialdow, p) )
	aff_transmg!(int) = begin 
		Kg -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :G, flinks, glinks, hr, carestage)
		)
	end
	j_transmg = VariableRateJump(rate_transmg, aff_transmg!; lrate=lrate_transmg, urate=hrate_transmg, rateinterval=rint)

	rate_transmh(u,p,t) = hr*transmissionrate(carestage, infstage, :H, t, tinf, tinfectious, initialdow, p) 
	hrate_transmh(u,p,t) = max(1.2*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p) )
	lrate_transmh(u,p,t) = min(0.8*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinf, tinfectious, initialdow, p) )
	aff_transmh!(int) = begin 
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t, tinf,initialdow), region, :H, flinks, glinks, hr, carestage)
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

	## hospital TODO add hospital -> discharged rate 
	ratehospital(u,p,t) = (( (carestage in (:undiagnosed,:GP)) & (severity in (:severe,)) )) ? p.hospadmitrate : 0.0 
	aff_hosp!(int) = begin 
		carestage = :admittedhospital  
		thospital = int.t 
	end 
	j_hospital = ConstantRateJump( ratehospital, aff_hosp! )

	## icu TODO add icu -> hosp rate  
	rateicu(u,p,t) = (( (carestage in (:undiagnosed,:GP,:admittedhospital)) & (severity in (:severe,)) )) ? p.icurate : 0.0 
	aff_icu!(int) = begin 
		carestage = :icu  
		ticu = int.t 
		sampled = (rand() < p.psampled )
	end 
	j_icu = ConstantRateJump( rateicu, aff_icu! )


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
		# @show homeregion 
		# @show commuteregion 
		# @show region
		# @show setdiff( [homeregion, commuteregion], [region] ) 
		# (homeregion!=commuteregion) && (region= setdiff( [homeregion, commuteregion], [region] )[1] )
		if (homeregion!=commuteregion) & ( region in [homeregion, commuteregion] )
			region= setdiff( [homeregion, commuteregion], [region] )[1]
		end
	end 
	j_commute = ConstantRateJump( ratecommute, aff_commute! )

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
	jumps = [ j_gainf
		, j_gaing
		, j_losef
		, j_loseg
		, j_transmf
		, j_transmg
		, j_transmh 
		, j_gp 
		, j_hospital 
		, j_icu
		# , j_migration
		, j_commute
	]
	# jdep complete graph works, but is probably slower 
	jdep = repeat([collect(1:length(jumps))],length(jumps)) # complete graph 
	simgenprob1 = JumpProblem(simgenprob0, Coevolve(), jumps...; dep_graph=jdep)
	sol = solve(simgenprob1, SSAStepper())

	# # check if sampled and generate tseqdeposition
	# sampled && (tsampled = rand(Uniform(ticu, tfin)) ) # TODO should have icu discharge distinct from recovery
	# sampled && ( tseqdeposition = tsampled + rand( Uniform( p.lagseqdblb , p.lagsseqdbub)) )
	
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
		, sol
		, tinf
		, tgp
		, thospital
		, ticu 
		# , tsampled 
		# , tseqdeposition 
		, trecovered 
		, contacttype
		, (flinks,glinks,hr)
		, iscommuter 
		, initialregion
		, homeregion
		, commuteregion
		, initialdow # day 1-7 when infection ocurred 
		, isnothing(donor) ? 0 : donor.generation+1 # generation 
	)
end

Infection(p, h::DataFrameRow, donor::Infection) = Infection(p; pid=h["pid2"]
							    , region=h["region"]
							, tinf=h["timetransmission"]
							, initialdow = h["dayofweek"]
							, contacttype=h["transmissiontype"]
							, donor=donor
)

# function Infection(p; pid = "0", region="TLI3", tinf = 0.0, initialdow = 1, contacttype = :nothing, donor::Union{Nothing,Infection} = nothing)

function simgeneration(p, prevgen::Array{Infection}; maxtime = Inf)
	length(prevgen) == 0 && return Array{Infection}([]) 
	Array{Infection}(
		[Infection(p, h, u) for u in prevgen for h in eachrow(u.H) if h["timetransmission"] < maxtime]
	)
end

function simtree(p; region="TLI3", initialtime=0.0, maxtime=30.0, maxgenerations::Int64=10, initialcontact=:H)
	# clustthreshold::Float64 = 0.005,
	
	@assert maxgenerations > 0
	@assert maxtime > initialtime

	simid = UUIDs.uuid1() |> string 

	g = Array{Infection}([ 
 	 	Infection(p; pid = "$(simid)-0", region=region, tinf = initialtime, contacttype= initialcontact) 
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
	
	dfargs = [(u.dpid, u.pid, u.tinf, u.contacttype, u.initialregion) for u in G if !ismissing(u.dpid)] # TODO initial region? 
	D = length(dfargs) > 0 ? 
		DataFrame(dfargs, [:donor, :recipient, :timetransmission, :contacttype, :region]) :  
		DataFrame([:donor => nothing, :recipient => nothing, :timetransmission => nothing, :contacttype => nothing, :region => nothing])
		
	dfargs1 = [(u.pid, u.tinf, u.tgp, u.thospital, u.ticu, u.trecovered, u.iscommuter, u.homeregion, u.commuteregion, u.generation,
				u.degree...) for u in G]
	Gdf = DataFrame(dfargs1
		, [:pid, :tinf, :tgp, :thospital, :ticu, :trecovered, :iscommuter, :homeregion, :commuteregion, :generation, :F, :G, :H ]
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
	# TODO should weight by pop density 
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
		simtree(p; region=r.region, initialtime=t, maxtime=maxtime, maxgenerations=maxgenerations, initialcontact=:H)
	end
	G = reduce( vcat, map(x-> x.G,trs))
	D = reduce( vcat, map(x-> x.D,trs))
	(; G = G, D = D, nimports = nimports, timports = timports  )
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

sampleforest(fo, p) = begin 
	# icu cases 
	G = fo.G[ isfinite.(fo.G.ticu), : ]
	# sample size 
	n = rand( Binomial( size(G,1), p.psampled ))
	# subforest 
	G1 = G[sample( 1:size(G,1), n, replace=false ), :]

	tsample = map( g -> rand(Uniform(g.ticu,g.trecovered)) , eachrow(G1) )
	G1.tsample = tsample 

	(; 
		tsample = sort( tsample  )
		, G = G1 
		, n = n 
		, firstsample = minimum(tsample)
	)
end
