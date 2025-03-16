#= 
# DONE 
# - change units to years (from days)
# - switch termination condition to time 
#	- make explicit start time end time 
# - make transmission rate func(t)  ; p.T scale should be dynamic 
# - introduce sequencing rate func(t)
# - make glinks parametric; tunable variance 
# - make oo links parameterised mean,variance
#
# TODO 
#
# ? inflate transmission in aids
# ? allow loss from care
# ? make diagnosis rate func(t) 
#
# Calibration
# - 1) sequencing 
#	* fix to given probability that a case is sequences over time 
# - diagnosis rate fixed?? 
# - 2) transmrate(t) ; Interpolations.jl spline? 
#	* match _growth_ rate in cases 
# - 3) glink and oo variance; trans prob scale; ν half-life of transm scale 
#	* match clustersize distribution 
# 2 & 3 both depend on variance of degree distribution, so will need to be calibrated simultaneously 
=#

using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using DataFrames
import UUIDs 
import StatsBase 
using Interpolations

import SpecialFunctions as SF 

using Plots 
default(;lw = 4)

using Debugger 

const UNDIAGNOSED = 0 
const DIAGNOSED = 1
const SUPPRESSED = 2 
const UNDIAGNOSED_EHI = 3 
const UNDIAGNOSED_CHRONIC = 4 
const UNDIAGNOSED_AIDS = 5 
const DECEASED = 6

const  MAXDURATION = 20.0

struct Infection
	pid::String # infection id 
	dpid::Union{String,Missing}
	H::DataFrame # transmission linelist 
	R::Int64 # number transmissions 
	sol::ODESolution
	tinf::Float64 # time infected 
	tdiagnosed::Float64 # time diagnosed 
	tsequenced::Float64 # time sequenced
	d::Float64 # distance to donor 
	contacttype::Symbol # cause of infection
	degree::Tuple{Int64,Int64,Float64}
	generation::Int64
end


# default parameters 
P = ( τ₀ = 0.001171 #  -log( 1 - .0075)* (7/1.5)*(1/30) ; O4Y4epigen  # initial transm prob per act (acute)
	, τ₁ = 0.0002646695 # O4Y4epigen # final transm prob per act {width=70%}
	# , T = 4.299657975412291 # scale for transmission prob 
	, ftransmscale =  interpolate([ 1990.0, 2010.0, 2020.0]
			       , [ 4.299657975412291, 4.299657975412291 , 4.299657975412291  ] 
			       , SteffenMonotonicInterpolation())
	, ν = -log(0.5)/(0.5) # decay rate transmission probability # half-life 6 months 
	, fcont = 1.5/7.0*365.0 # contact rate for flinks  
	, gcont = 1.0/7.0*365.0 # contact rate for glinks 
	, δ = -365*log( 1.0 - 0.0191 ) / 30.0 # ~1/4yrs # diagnosis rate; O4Y4epigen *
	, κ = 365*0.004735967  # ~1/211 days; ( (-log( 1 - 0.2367 ) / 30 )^(-1)+(-log(1-0.2590)/30)^(-1) )^(-1) O4Y4epigen , combining care and suppression intervals 
	, γ_ehi = 1.0 # volz pmed 2013 
	, γ_chron = (2.1) # volz pmed 2013 
	, shape_chron = 3.0
	, γ_aids = 1.0/(2.0)
	# , α = -365*log(1.0-781.0/1049.0)/30.0 # sequencing rate following diagnosis; O4Y4epigen/sd_sequence_prob_report 
	# , fα = interpolate([1990.0, 2010.0, 2020.0]
	# 		, [16.6,16.6,16.6]
	# 	     , SteffenMonotonicInterpolation() ) 
	, fpsequenced= interpolate([1990.0, 2010.0, 2020.0]
		, [.5, .5, .5]
		, SteffenMonotonicInterpolation() ) 
	, ρ_diagnosed = 0.50 # TODO # transmission reduction 
	, ρ_suppressed = 0.0 # transmission reduction 
	, frate = 365*1.0/7.0/30.0 # rate of gaining & losing flinks
	, grate = 365*1.0/2.0/30.0 # rate of gaining and losing
	, fdist = [0.57279772, 0.40783034, 0.01937194]
	, fdist1 = [0.57279772, 0.40783034, 0.01937194].*collect(0:2) / sum([0.57279772, 0.40783034, 0.01937194].*collect(0:2) ) # excess degree 
	, gdist =  [0.60012235,0.20575041,0.09726754,0.09685971]  
	, gdist1 = [0.60012235,0.20575041,0.09726754,0.09685971] .* collect(0:3) / sum([0.60012235,0.20575041,0.09726754,0.09685971] .* collect(0:3))  # excess degree
	, gnegbinomr = 1.5995135923340424 # based on gdist  
	, gnegbinomp = 0.6983622109747697
	, oorateshape = 0.26298  
	, ooratescale = 0.07605*365
	, oorateshape1 = 1.26298 # excess rate distribution (also gamma)
	, ooratescale1 = 0.07605*365
	, varinflation = 1 # inflates variance by this factor
	, μ = 0.001 # mean clock rate -- additive relaxed clock
	, ω = 0.5 # variance inflation
)


function sampdegree(p; contacttype = :nothing)
	gnegbinomr = p.gnegbinomr*p.gnegbinomp/p.varinflation
	gnegbinomp = p.gnegbinomp/(p.gnegbinomp + p.varinflation*(1-g.negbinomp))
	r = Gamma( p.oorateshape/p.varinflation, p.ooratescale*p.varinflation ) |> Base.rand 
	kf = StatsBase.wsample( 0:2, p.fdist  )
	kg = rand( NegativeBinomial(gnegbinomr, gnegbinomp) )
	if contacttype == :F # note this counts the link that transmitted infection
		kf = StatsBase.wsample(0:2, p.fdist1 )
	elseif contacttype == :G
		# kg = StatsBase.wsample( 0:3, p.gdist1 )
		kg = rand( NegativeBinomial(gnegbinomr+1, gnegbinomp) ) + 1 
	elseif contacttype == :H 
		r = Gamma( p.oorateshape1/p.varinflation, p.ooratescale1*p.varinflation ) |> Base.rand 
	end 
	[ kf, kg, r ]	
end

function transmissionprob(carestage, t,tinf,p ) 
	ρ = 1.0 
	if carestage == DIAGNOSED
		ρ = p.ρ_diagnosed
	elseif carestage == SUPPRESSED 
		ρ = p.ρ_suppressed
	elseif carestage == DECEASED
		ρ = 0.0 
	end 
	o = p.τ₁ + ( p.τ₀-p.τ₁ ) * exp( -p.ν * t )
	# TODO dynamic transm here, tinf,  
	transmscale = p.ftransmscale( tinf + t )
	# p.T * ρ * o # scale * suppression * f(vl) 
	transmscale * ρ * o # scale * suppression * f(vl) 
end

function simgendist( t0, t1, p ; s = 1000)
	@assert t1 >= t0 
	muts = Base.rand( 
		NegativeBinomial(s * p.μ * (t1-t0) / p.ω,  1.0 / ( 1 + p.ω ) ) # NOTE prob = 1-<paramter in ARC paper>
		)
	float( muts ) / s 
end

function simgendist( tuv, p; s = 1000 )
	# Gamma(  p.μ * (tuv)/(1.0+p.ω),  1.0+p.ω ) |> Base.rand
	simgendist( 0., tuv, p; s = s)
end

function Infection(p ; pid = "0", tinf = 0.0, contacttype = :nothing
		   , donor::Union{Nothing,Infection} = nothing
	)
	H = DataFrame(pid1 = String[], pid2 = String[]
	       , timetransmission = Float64[]
	       , transmissiontype = Symbol[]
	       # , gendist0 = Float64[]
	       , degreef = Int64[]
	       , degreeg = Int64[] 
	       , oorate = Float64[] 
	)
	flinks, glinks, hr = sampdegree(p; contacttype = contacttype )
	u₀ = [ Poisson( flinks ) |> Base.rand 
		, Poisson( glinks ) |> Base.rand 
		, UNDIAGNOSED_EHI 
	]
	if contacttype == :F 
		u₀[1] =  max(0, u₀[1]-1 )
	elseif contacttype == :G
		u₀[2]  = max(0, u₀[2]-1 )
	end 
	tspan = (tinf, tinf + MAXDURATION)
	tseq = Inf # time of sequencing 
	tdiagnosed = tinf + MAXDURATION # time of diagnosis 
	R = 0 #cumulative transm 
	diagnosedwithaids = false 

	Γchron = Gamma( p.shape_chron,  p.γ_chron )
	gamhazard( t) =	pdf( Γchron, t ) / ( 1.0-cdf(Γchron, t ))

	kv = sampdegree(p)
	
	rate_gainf(u, p, t) = flinks * p.frate
	rate_losef(u, p, t) = u[1] * p.frate 
	
	rate_gaing(u, p, t) = glinks * p.grate
	rate_loseg(u, p, t) = u[2] * p.grate 
	
	rate_ehi2chron(u,p,t) = u[3]==UNDIAGNOSED_EHI ?  p.γ_ehi : 0.0
	tchronstart = missing
	rate_chron2aids(u,p,t) = u[3]==UNDIAGNOSED_CHRONIC ? gamhazard( t-tchronstart ) : 0.0 
		hrate_chron2aids(u,p,t) = 1
		lrate_chron2aids(u,p,t) = 0.0
	rate_aids2death(u,p,t) = u[3]==UNDIAGNOSED_AIDS ? p.γ_aids : 0.0 
	rate_diagnosis(u, p, t) = u[3]∈[UNDIAGNOSED_EHI,UNDIAGNOSED_CHRONIC,UNDIAGNOSED_AIDS] ?  p.δ*(1+4*(u[3]==UNDIAGNOSED_AIDS)) : 0.0 
	rate_care(u,p,t) = u[3]==DIAGNOSED ? p.κ*(1+4diagnosedwithaids) : 0.0  

	rate_transmf(u,p,t) = u[1] * p.fcont * transmissionprob(u[3], t-tinf, tinf, p )
	hrate_transmf(u,p,t) = u[1] * p.fcont * p.τ₀ * p.ftransmscale(t)  #p.T
		lrate_transmf(u,p,t) = 0.0 # u[1] * p.fcont * p.τ₁
	rate_transmg(u,p,t) = u[2] * p.gcont * transmissionprob(u[3], t-tinf, p )
		hrate_transmg(u,p,t) = u[2] * p.gcont * p.τ₀ * p.ftransmscale(t)   
		lrate_transmg(u,p,t) = 0. # u[2] * p.gcont * p.τ₁
	rate_transmh(u,p,t) = hr * transmissionprob(u[3], t-tinf,p)
		hrate_transmh(u,p,t) = hr * p.τ₀ * p.ftransmscale(t) 
		lrate_transmh(u,p,t) = 0.0 # hr * p.τ₁
	rint(u,p,t) = 1.0 #Inf 

	function aff_gainf!(int)
		int.u[1] += 1 
	end

	function aff_gaing!(int)
		int.u[2] += 1
	end

	function aff_losef!(int)
		int.u[1] -= 1
	end

	function aff_loseg!(int)
		int.u[2] -= 1
	end

	function aff_ehi2chron!(int)
		int.u[3] = UNDIAGNOSED_CHRONIC 
		tchronstart = int.t 
	end 

	function aff_chron2aids!(int)
		int.u[3] = UNDIAGNOSED_AIDS 
	end

	function aff_aids2death!(int)
		tdiagnosed = int.t
		int.u[3] = DECEASED
	end

	function aff_diagnosis!(int)
		# α = p.fα( int.t ) 
		α = p.α
		pseq = p.fpsequenced(int.t)
		tseq = Inf
		if rand() < pseq:
			tseq = int.t + (Exponential(1.0/α) |> Base.rand  )
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

	function aff_transmf!(int)
		int.u[1] -= 1
		R +=  1
		nextpid = pid * ".$(R)"
		push!( H
			, (pid, nextpid, int.t, :F
			# , simgendist(tinf, int.t, p)
			, flinks, glinks, hr ) 
		)
	end
	
	function aff_transmg!(int)
		int.u[2] -= 1
		R +=  1
		nextpid = pid * ".$(R)"
		push!( H
			, (pid, nextpid, int.t, :G 
			# , simgendist(tinf, int.t, p)
			, flinks, glinks, hr ) 
		)
# @show (R, int.t,int.u,  )
	end
	
	function aff_transmh!(int)
		R +=  1
		nextpid = pid * ".$(R)"
		push!( H
			, (pid, nextpid, int.t, :H
			# , simgendist(tinf, int.t, p)
			, flinks, glinks, hr ) 
		)
	end
	
	j_gainf = ConstantRateJump( rate_gainf, aff_gainf!) 
	j_gaing = ConstantRateJump( rate_gaing, aff_gaing!) 

	j_losef = ConstantRateJump( rate_losef, aff_losef!) 
	j_loseg = ConstantRateJump( rate_loseg, aff_loseg!) 

	j_ehi2chron = ConstantRateJump(rate_ehi2chron, aff_ehi2chron! )
	j_aids2death = ConstantRateJump(rate_aids2death, aff_aids2death! )
	j_diagnosis = ConstantRateJump( rate_diagnosis, aff_diagnosis! )
	j_care = ConstantRateJump( rate_care, aff_care! ) 

	j_chron2aids = VariableRateJump( rate_chron2aids, aff_chron2aids!; lrate= lrate_chron2aids, urate =hrate_chron2aids, rateinterval= rint) # 
	j_transmf = VariableRateJump( rate_transmf, aff_transmf!; lrate= lrate_transmf, urate =hrate_transmf, rateinterval= rint) # 
	j_transmg = VariableRateJump( rate_transmg, aff_transmg!; lrate = lrate_transmg, urate =hrate_transmg, rateinterval= rint )
	j_transmh = VariableRateJump( rate_transmh, aff_transmh!; lrate = lrate_transmh, urate = hrate_transmh, rateinterval = rint )
	
	simgenprob0 = DiscreteProblem( u₀, tspan, p )
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

	# 
	# jdep = [ # when i'th jump occurs, these jump rates need to be recalculated 
	# 	[3, 10 ] #					1 gainf 
	# 	, [4, 11 ] #					2 gaing 
	# 	, [3, 10 ] #					3 losef 
	# 	, [4, 11 ] #					4 loseg 
	# 	, [5, 6] #					5 ehi2chron 
	# 	, [6, 7] #					6 chron2aids
	# 	, [7,8,10,11,12] #				7 aids2death 
	# 	, [5,6,7,8,9,10,11,12] #			8 diagnosis 
	# 	, [8,9,10,11,12] #				9 care 
	# 	, [3,10] #					10 transmf
	# 	, [4,11] #					11 transmg 
	# 	, [] #						12 transmh 
	# ]
	
	# jdep complete graph works, but is probably slower 
	jdep = repeat( [collect( 1:12 )] , 12 ) # complete graph 

	simgenprob1 = JumpProblem( simgenprob0, Coevolve(), jumps...; dep_graph = jdep )
	
	sol = solve( simgenprob1, SSAStepper() )

	if isnothing(donor) || isinf( tseq )
		d = Inf 
	elseif isinf( donor.tsequenced )
		d = Inf 
	else 
		d = simgendist( abs((tseq-tinf) + (donor.tsequenced-tinf)), p  )
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
     	);
end

Infection(p, h::DataFrameRow, donor::Infection ) = Infection(p; pid = h["pid2"], tinf = h["timetransmission"], contacttype=h["transmissiontype"], donor = donor );


function simgeneration(p, prevgen::Array{Infection})
	length( prevgen ) == 0 && return Array{Infection}([]) 
	Array{Infection}(
		[ Infection(p, h, u) for u in prevgen for h in eachrow(u.H)  ]
	)
end

function simbp(p ; initialtime=1990.0,  maxtime = 2020.0, maxgenerations::Int64 = 100, initialcontact=:G)
	# clustthreshold::Float64 = 0.005,
	
	@assert maxgenerations > 0

	g = Array{Infection}( [ 
 	 	 Infection(p; pid = "0", tinf = initialtime, contacttype= initialcontact );
	] )

	G = g 
	H = g[1].H 
	for igen in 2:maxgenerations
		g = simgeneration(p, g )
		g = [ infection for infection in g if infection.tinf < maxtime ]
		if  length( g ) > 0 
			H = vcat( H, vcat( [x.H for x in g]... ) )
			try
				push!( G, g...)
			catch
				@bp 
			end
		end

	end
	dfargs = [ (u.dpid, u.pid, u.d, u.tinf, u.contacttype) for u in G if !ismissing(u.dpid) ]
	D = length(dfargs)>0 ? 
		DataFrame( dfargs,  [:donor, :recipient, :distance, :timetransmission, :contacttype]  ) :  
		DataFrame( [:donor => nothing, :recipient => nothing, :distance => nothing, :timetransmission => nothing, :contacttype => nothing] )
		
	dfargs1 = [ (u.pid, u.tsequenced, u.tdiagnosed, u.tinf, u.generation
	, u.degree...) for u in G ]
	Gdf = DataFrame( dfargs1, [:pid, :timesequenced, :timediagnosed, :timeinfected, :generation, :Fdegree, :Gdegree, :Hdegree ] )
	
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

