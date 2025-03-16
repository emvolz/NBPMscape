using NBPMscape
using Test

@testset "NBPMscape.jl" begin
	# Test infection creation
	infection = Infection(NBPMscape.P; pid="test", tinf=2000.0, contacttype=:G)
	@test infection.pid == "test"
	@test infection.tinf == 2000.0
	@test infection.contacttype == :G
	
	# Test simulation
	result = simbp(NBPMscape.P; initialtime=2000.0, maxtime=2001.0, maxgenerations=2)
	@test result.G isa DataFrame
	@test result.D isa DataFrame
	@test !isempty(result.infections)
end
