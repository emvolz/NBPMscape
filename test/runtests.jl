using NBPMscape
using Test

if false 

	using NBPMscape
	using Debugger
	using Revise
	using Test
	using StatsBase
	using Distributions

	tr = simtree(NBPMscape.P; initialtime=1.0, maxtime=50.0, maxgenerations=100)

	# ~ 2.25 sec
	fo = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100) ;

	fo = simforest(NBPMscape.P; initialtime=0.0, maxtime=120.0, maxgenerations=100); [ fo.nimports , size(fo.G ,1)]

	P = merge( NBPMscape.P, (;nimports = 1000))
	fo = simforest(P; initialtime=0.0, maxtime=60.0, maxgenerations=100); [ fo.nimports , size(fo.G ,1)]

	G = fo.G[ isfinite.(fo.G.ticu), : ]
	psample = .05 
	n = rand( Binomial( size(G,1), psample ))
	G1 = G[sample( 1:size(G,1), n, replace=false ), :]

	
end 


if false 

	using NBPMscape
	using Debugger
	using Revise
	using Test
	using StatsBase
	using Distributions

	d = NBPMscape.CAAIMPORTS
	region = wsample( d.ITL225CD, d.pax_2024_per_day )
	prd = deepcopy( NBPMscape.COMMUTEINPROB[region] )
	("na" in prd.index2name) && (delete!( prd, "na" ))
	wsample( prd.index2name, prd.data )

	sampleimportregion( P)

end 


if false 

	using NBPMscape
	using Debugger
	using Revise
	infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)

	infectivitytoR( 5.0, nsims = 10 ) isa Union{Missing,Real}

	result = simtree(NBPMscape.P; initialtime=1.0, maxtime=50.0, maxgenerations=2)


	ch = simtree(NBPMscape.P; initialtime=1.0, maxtime=50.0, maxgenerations=100)

end 

if false
	using Revise
	using Debugger
	using StatsBase

	# @run infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)
	infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)

	P = NBPMscape.P; 
	P = merge(P, (;infectivity=128.0));
	transmissionrate(:undiagnosed, :infectious, :G, 3.0, 2.0, 1, P) 

	P = NBPMscape.P
	infs = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:G), 1:1000);
	H = reduce(vcat, [x.H for x in infs] )
	H.timetransmission |> summarystats
	[ x.R for x in infs ] |> summarystats

	reduce( .+ , map( i->sampdegree(P;contacttype = :G), 1:1000)) ./ 1e3
	reduce( .+ , map( i->sampdegree(P;contacttype = :F), 1:1000)) ./ 1e3

	nsims = 1000
	@time foreach(_->Infection(P; pid="test", tinf=0.0, contacttype=:F), 1:nsims);
  	# 0.459718 seconds (4.29 M allocations: 340.712 MiB, 14.09% gc time, 2.18% compilation
	# 0.00046
	
	nsims = 1000
	infsF = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:F), 1:nsims);
	infsG = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:G), 1:nsims);
	infsH = map(i->Infection(P; pid="$(i)", tinf=0.0, contacttype=:H), 1:nsims);
	inftoR(inf) = begin 
		cmh = countmap( inf.H.transmissiontype )
		map( k-> (k in keys(cmh)) ? cmh[k] : 0, [:F, :G, :H] )
	end
	infstoR(infs) = reduce( .+, map(inftoR,infs));
	RF = infstoR(infsF) ./ nsims ;
	RG = infstoR(infsG) ./ nsims ;
	RH = infstoR(infsH) ./ nsims ;
	NGM = [ RF RG RH ]

	using LinearAlgebra
	eigvals(NGM)
	eigvecs(NGM)

end




@testset "NBPMscape.jl" begin
	# Test infection creation
	infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)
	@test infection.pid == "test"
	@test infection.tinf == 1.2
	@test infection.contacttype == :G
	
	@test  infectivitytoR( 5.0, nsims = 10 ) isa Union{Missing,Real}

	# Test simulation of single seed 
	result = simtree(NBPMscape.P; initialtime=1.0, maxtime=50.0, maxgenerations=2)
	@test result.G isa DataFrame
	@test result.D isa DataFrame
	@test !isempty(result.infections)

	# test multiple imports 
	forest = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=2); 
	@test forest.G isa DataFrame
	@test forest.D isa DataFrame
	@test forest.nimports > 0
end

