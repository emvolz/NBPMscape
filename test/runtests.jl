using NBPMscape
using Test

if false
	using Revise
	using Debugger
	using StatsBase

	# @run infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)
	infection = Infection(NBPMscape.P; pid="test", tinf=1.2, contacttype=:G)

	P = NBPMscape.P; 
	P = merge(P, (;infectivity=128.0));
	transmissionrate(:undiagnosed, :infectious, :G, 3.0, 2.0, 1, P) 

	infs = map(i->Infection(P; pid="$(i)", tinf=1.2, contacttype=:G), 1:100);
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
	
	# # Test simulation
	# result = simbp(NBPMscape.P; initialtime=2000.0, maxtime=2001.0, maxgenerations=2)
	# @test result.G isa DataFrame
	# @test result.D isa DataFrame
	# @test !isempty(result.infections)

end

