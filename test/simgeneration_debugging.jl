### Debugging of simgeneration function
#=
Previously there was an issue with the simgeneration function in core.jl.
However, following the incorporation of age disaggregation of severity, 
the error that was being generated is not being generated as at 16 Sep 2025.
The scripts below contain the comments previously made and an Alternative
version of the function used for debugging.

=#

function simgeneration(p, prevgen::Array{Infection}; maxtime = Inf)
	length(prevgen) == 0 && return Array{Infection}([])
	# TODO THERE IS AN ERROR HERE - possibly in "if h["timetransmission"] < maxtime" which apparently sometimes
	# tries to compare a vector with a single value
	# Loops through the each Infection u in the prevgen array of Infections and loops through each row h in u.H.
	# If the timetransmission value in row h in H in Infection u is less than maxtime then a new Infection object 
	# will be created using parameter p (the parameter set), h (row in H) and u (the prevgen Infection) 
	Array{Infection}(
		[Infection(p, h, u) for u in prevgen for h in eachrow(u.H) if h["timetransmission"] < maxtime]
	)
end

#TEST
# Alternative version of simgeneration() for debugging
#prevgen = g
#maxtime=60.0
#p=P
function simgeneration_alt(p, prevgen::Array{Infection}; maxtime = Inf)
	#Test
	#println("length of g = ",length(g))
	#println(prevgen[1])
	length(prevgen) == 0 && return Array{Infection}([])	
	newgen = Infection[]
	for u in 1:length(prevgen)
		for h in 1:nrow(prevgen[u].H)
			#println("u=",u", ","h=",h)
			println("Infection ",u," H row",h," has time of transmission = ",prevgen[u].H[h,:]["timetransmission"])
			println("t_transmission",prevgen[u].H[h,:]["timetransmission"])
			if prevgen[u].H[h,:]["timetransmission"] < maxtime
				#Array{Infection}
				#push!( Array{Infection},  Infection(p,h,u) )
				push!( newgen,  Infection(p,prevgen[u].H[h,:],prevgen[u]) )
			end
		end
	end
	return( newgen )
end