default(;lw = 4)

const MAXDURATION = 180.0

struct Infection
	pid::String # infection id 
	dpid::Union{String,Missing}
	H::DataFrame # transmission linelist 
	R::Int64 # number transmissions 
	sol::ODESolution
	tinf::Float64 # time infected 
	tgp::Float64 # time 
	thospital::Float64 # time 
	ticu::Float64 # time 
	tsampled::Float64 # time
	tseqdeposition::Float64 # time for sequencing + qc + bioninf 
	trecovered::Float64 
	d::Float64 # distance to donor 
	contacttype::Symbol # cause of infection
	degree::Tuple{Int64,Int64,Float64}
	initialregion::Int
	finalregion::Int
	initialdayofweek::Int # day 1-7 when infection ocurred 
	generation::Int64
end

const AGEGROUPS = (:youth, :youngadult, :adult, :old, :elderly) 
const CARE  = (:undiagnosed, :GP, :admittedhospital, :admittedicu, :discharged) 
const SEVERITY = (:mildorasymptomatic, :moderate, :severe )
const STAGE = (:latent, :infectious, :recovered)

const CONTACTTYPES = (:F, :G, :H)

const NREGIONS = 50
const REGIONS = ( :TODO )

# TODO temporary migration matrix 
MIGMATRIX = fill( 1.0, NREGIONS, NREGIONS )
for i in 1:NREGIONS
	MIGMATRIX[i,i] = 100.0 
	MIGMATRIX[i,:] ./= sum( MIGMATRIX[i,:] )
end


# default parameters 
P = ( 
	#τ₀ = 0.1 #  -log( 1 - .0075)* (7/1.5)*(1/30) ; # initial transm prob per act (acute)
	#, τ₁ = 0.01 # # final transm prob per act 
	# , T = 4.299657975412291 # scale for transmission prob 
	# , ftransmscale =  interpolate([ 1990.0, 2010.0, 2020.0]
	# 			   , [ 4.299657975412291, 4.299657975412291 , 4.299657975412291  ] 
	# 			   , SteffenMonotonicInterpolation())
	# , ν = -log(0.5)/(0.5) # decay rate transmission probability # half-life 6 months 
	
	fcont = 1.0 # relative contact rate for flinks  (household) # TODO get from contact tracing studies 
	, gcont = .25 # contact rate for glinks (work)
	, oocont = .05 # contact rate stranger 

	, dowcont = (.2, .75, 1.0, 1.0, 1.0, .8, .3) # day of week scale Sun-Sat 

	, infectivity = 1.0 # scales transmission rate 
	, infectivity_shape = 2.2 # Metcalf Nat Comm 2021 
	, infectivity_scale = 2.5 

	, latent_shape = 2.3 # Casey 2021 ?
	, latent_scale = 1.2
	, infectious_shape = 2.1 # Verity 2020 , mean 22 d ?
	, infectious_scale = 10.57 

	# , δ = -365*log( 1.0 - 0.0191 ) / 30.0 # ~1/4yrs # diagnosis rate; O4Y4epigen *
	# , κ = 365*0.004735967  # ~1/211 days; ( (-log( 1 - 0.2367 ) / 30 )^(-1)+(-log(1-0.2590)/30)^(-1) )^(-1) O4Y4epigen , combining care and suppression intervals 
	# , γ_inf = 2.1*365 # 
	# , shape_inf = 3.0*365 # 
	
	, ρ = 0.250 #  transmission reduction 
	
	, frate = 0.0 # rate of gaining & losing flinks
	, grate = 1/30.0 # rate of gaining and losing

	, fnegbinomr = 8.1 # Approximate household size distribution TODO 
	, fnegbinomp = 0.86 
	, gnegbinomr =  3 # Approximate workplace size distribution TODO 
	, gnegbinomp = 0.25 
	, oorateshape = 1.975 #   # Approximate stranger contacts (e.g. public transport) TODO 
	, ooratescale = 5.0 # 
	, oorateshape1 = 2.9750  # excess rate distribution (also gamma)
	, ooratescale1 = 5.0 

	, μ = 0.001 # mean clock rate -- additive relaxed clock
	, ω = 0.5 # variance inflation
)


function sampdegree(p; contacttype = :nothing)
	gnegbinomr = p.gnegbinomr*p.gnegbinomp/p.varinflation
	gnegbinomp = p.gnegbinomp/(p.gnegbinomp + p.varinflation*(1-g.negbinomp))
	r = Gamma( p.oorateshape/p.varinflation, p.ooratescale*p.varinflation ) |> Base.rand 
	kf = rand( NegativeBinomial(fnegbinomr, fnegbinomp) )
	kg = rand( NegativeBinomial(gnegbinomr, gnegbinomp) )
	if contacttype == :F # note this counts the link that transmitted infection
		kf = rand( NegativeBinomial(fnegbinomr+1, fnegbinomp) ) + 1 
	elseif contacttype == :G
		kg = rand( NegativeBinomial(gnegbinomr+1, gnegbinomp) ) + 1 
	elseif contacttype == :H 
		r = Gamma( p.oorateshape1/p.varinflation, p.ooratescale1*p.varinflation ) |> Base.rand 
	end 
	[ kf, kg, r ]	
end

function dayofweek(t, initialdow)
	Int( floor( t-initialdow-1 ) % 7  ) + 1
end

function transmissionrate(carestage, infstage, contacttype, t, tinfectious, initialdow, p) 
	dow = dayofweek(t,initialdow) 
	ρ = 1.0 
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

function Infection(p; pid = "0", tinf = 0.0, initialdow = 1, contacttype = :nothing, donor::Union{Nothing,Infection} = nothing)
	# transmission line list 
	H = DataFrame(pid1 = String[], pid2 = String[]
	       , timetransmission = Float64[]
	       , transmissiontype = Symbol[]
	       # , gendist0 = Float64[]
	       , degreef = Int64[]
	       , degreeg = Int64[] 
	       , oorate = Float64[] 
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
	carestage = :undiagnosed 
	infstage = :latent 
	agegroup = :adult # TODO sample 
	severity = :moderate # TODO sample 
	region = :london # TODO sample 

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
	tspan = (tinfectious, tfin)
	

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
	rate_transmf(u,p,t) = Kf*transmissionrate(carestage, infstage, :F, t, tinfectious, initialdow, p) 
	hrate_transmf(u,p,t) = max(1.2*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3, tinfectious, initialdow, p) )
	lrate_transmf(u,p,t) = min(0.8*rate_transmf(u,p,t), Kf*transmissionrate(carestage, infstage, :F, t+3, tinfectious, initialdow, p) )
	aff_transmf!(int) = begin 
		Kf -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t,initialdow), :F, flinks, glinks, hr, carestage)
		)
	end
	j_transmf = VariableRateJump(rate_transmf, aff_transmf!; lrate=lrate_transmf, urate=hrate_transmf, rateinterval=rint) # 
	
	rate_transmg(u,p,t) = Kg*transmissionrate(carestage, infstage, :G, t, tinfectious, initialdow, p) 
	hrate_transmg(u,p,t) = max(1.2*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3, tinfectious, initialdow, p) )
	lrate_transmg(u,p,t) = min(0.8*rate_transmg(u,p,t), Kg*transmissionrate(carestage, infstage, :G, t+3, tinfectious, initialdow, p) )
	aff_transmg!(int) = begin 
		Kg -= 1
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t,initialdow), :G, flinks, glinks, hr, carestage)
		)
	end
	j_transmg = VariableRateJump(rate_transmg, aff_transmg!; lrate=lrate_transmg, urate=hrate_transmg, rateinterval=rint)

	rate_transmh(u,p,t) = hr*transmissionrate(carestage, infstage, :H, t, tinfectious, initialdow, p) 
	hrate_transmh(u,p,t) = max(1.2*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinfectious, initialdow, p) )
	lrate_transmh(u,p,t) = min(0.8*rate_transmh(u,p,t), hr*transmissionrate(carestage, infstage, :H, t+3, tinfectious, initialdow, p) )
	aff_transmh!(int) = begin 
		R += 1
		nextpid = pid * ".$(R)"
		push!(H
			, (pid, nextpid, int.t, dayofweek(int.t,initialdow), :H, flinks, glinks, hr, carestage)
		)
	end
	j_transmh = VariableRateJump(rate_transmh, aff_transmh!; lrate=lrate_transmh, urate=hrate_transmh, rateinterval=rint)


	# care/sampling pathways TODO 
	rate_ehi2chron(u,p,t) = u[3]==UNDIAGNOSED_EHI ? p.γ_ehi : 0.0
	tchronstart = missing
	rate_chron2aids(u,p,t) = u[3]==UNDIAGNOSED_CHRONIC ? gamhazard(t-tchronstart) : 0.0 
		hrate_chron2aids(u,p,t) = 1
		lrate_chron2aids(u,p,t) = 0.0
	rate_aids2death(u,p,t) = u[3]==UNDIAGNOSED_AIDS ? p.γ_aids : 0.0 
	rate_diagnosis(u, p, t) = u[3]∈[UNDIAGNOSED_EHI,UNDIAGNOSED_CHRONIC,UNDIAGNOSED_AIDS] ? p.δ*(1+4*(u[3]==UNDIAGNOSED_AIDS)) : 0.0 
	rate_care(u,p,t) = u[3]==DIAGNOSED ? p.κ*(1+4diagnosedwithaids) : 0.0  

	rint(u,p,t) = 1.0 #Inf 

	function aff_diagnosis!(int)
		# α = p.fα(int.t) 
		α = p.α
		pseq = p.fpsequenced(int.t)
		tseq = Inf
		if rand() < pseq
			tseq = int.t + (Exponential(1.0/α) |> Base.rand)
		end
		tdiagnosed = int.t 
		if int.u[3] == UNDIAGNOSED_AIDS
			# will influence subsequent care rate 
			diagnosedwithaids = true 
		end 
		int.u[3] = DIAGNOSED
	end

	function aff_care!(int)
		int.u[3] = SUPPRESSED
	end

	j_ehi2chron = ConstantRateJump(rate_ehi2chron, aff_ehi2chron!)
	j_aids2death = ConstantRateJump(rate_aids2death, aff_aids2death!)
	j_diagnosis = ConstantRateJump(rate_diagnosis, aff_diagnosis!)
	j_care = ConstantRateJump(rate_care, aff_care!) 

	j_chron2aids = VariableRateJump(rate_chron2aids, aff_chron2aids!; lrate= lrate_chron2aids, urate=hrate_chron2aids, rateinterval=rint) # 
	
	simgenprob0 = DiscreteProblem(u₀, tspan, p)
	jumps = [ j_gainf
		, j_gaing
		, j_losef
		, j_loseg
		, j_ehi2chron
		, j_chron2aids
		, j_aids2death 
		, j_diagnosis
		, j_care
		, j_transmf
		, j_transmg
		, j_transmh 
	]

	# jdep complete graph works, but is probably slower 
	jdep = repeat([collect(1:12)], 12) # complete graph 

	simgenprob1 = JumpProblem(simgenprob0, Coevolve(), jumps...; dep_graph=jdep)
	
	sol = solve(simgenprob1, SSAStepper())

	if isnothing(donor) || isinf(tseq)
		d = Inf 
	elseif isinf(donor.tsequenced)
		d = Inf 
	else 
		d = simgendist(abs((tseq-tinf) + (donor.tsequenced-tinf)), p)
	end

	Infection(
		pid
		, isnothing(donor) ? missing : donor.pid 
		, H 
		, R
		, sol
		, tinf
		, tdiagnosed 
		, tseq 
		, d 
		, contacttype
		, (flinks,glinks,hr)
		, isnothing(donor) ? 0 : donor.generation+1
	)
end

Infection(p, h::DataFrameRow, donor::Infection) = Infection(p; pid=h["pid2"], tinf=h["timetransmission"], contacttype=h["transmissiontype"], donor=donor)


function simgeneration(p, prevgen::Array{Infection})
	length(prevgen) == 0 && return Array{Infection}([]) 
	Array{Infection}(
		[Infection(p, h, u) for u in prevgen for h in eachrow(u.H)]
	)
end

function simbp(p; initialtime=1990.0, maxtime=2020.0, maxgenerations::Int64=100, initialcontact=:G)
	# clustthreshold::Float64 = 0.005,
	
	@assert maxgenerations > 0

	g = Array{Infection}([ 
 	 	 Infection(p; pid = "0", tinf = initialtime, contacttype= initialcontact)
	])

	G = g 
	H = g[1].H 
	for igen in 2:maxgenerations
		g = simgeneration(p, g)
		g = [infection for infection in g if infection.tinf < maxtime]
		if length(g) > 0 
			H = vcat(H, vcat([x.H for x in g]...))
			G = vcat(G, g)
		end
	end
	
	dfargs = [(u.dpid, u.pid, u.d, u.tinf, u.contacttype) for u in G if !ismissing(u.dpid)]
	D = length(dfargs) > 0 ? 
		DataFrame(dfargs, [:donor, :recipient, :distance, :timetransmission, :contacttype]) :  
		DataFrame([:donor => nothing, :recipient => nothing, :distance => nothing, :timetransmission => nothing, :contacttype => nothing])
		
	dfargs1 = [(u.pid, u.tsequenced, u.tdiagnosed, u.tinf, u.generation,
				u.degree...) for u in G]
	Gdf = DataFrame(dfargs1, [:pid, :timesequenced, :timediagnosed, :timeinfected, :generation, :Fdegree, :Gdegree, :Hdegree])
	
	simid = UUIDs.uuid1() |> string 
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
