using Revise
using NBPMscape 
using JLD2 
using StatsBase
using CairoMakie
using CMPlot 
using PlotlyJS
using DataFrames


@load "covidlike-1.0-sims.jld2" 

sampps = [ .05, .1, .15 ] 

# ribbon plot 
times = range( 0.0, 60.0, length = 200 ) 	
it = reduce(hcat,  map( sims ) do s 
	map( t-> sum( s.G.tinf .< t ), times )
end 
)
bdata = mapslices( i->quantile(i, [.5, .025, .975] ), it, dims=2) 

fband = Figure(resolution=(3.5*72, 3.25*72), fontsize=11);
ax = Axis(fband[1, 1], xlabel="Days after first introduction", ylabel="Cumulative infections", title=""
	, yscale=log10
	  , limits = (nothing,(minimum(bdata[bdata[:,2].>0,2])*0.9, maximum(bdata[:,3])*1.1))
	);
band!(ax, times, bdata[:,2], bdata[:,3], color=(:blue, .15));
lines!(ax, times, bdata[:,1], color=:black);
fband
save( "covidlike-1.0-fband.pdf", fband ) 



# time to detect  
using Debugger
sp = .1 
spstats = map(sims)  do s 
	# try 
		ss = sampleforest( s , sp )
		t0 = ss.firstsample
		t3 = (ss.n>2) ? ss.tsample[3] : missing
		i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
	# catch 
# @bp 
	# end
	[ t0, t3, i0 ]
end
spstats = reduce( hcat, spstats )

sp = .25
spstats25 = map(sims)  do s 
	ss = sampleforest( s , sp )
	t0 = ss.firstsample
	t3 = (ss.n>2) ? ss.tsample[3] : missing
	i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
	[ t0, t3, i0 ]
end
spstats25 = reduce( hcat, spstats25 )

spdf = DataFrame( Days = vec(spstats')
	  	 , Count = repeat( [Symbol("1 case"), Symbol("3 cases"), :Infections], inner=size(spstats,2)))
spdf1 = spdf[ spdf.Count .!= :Infections, :]
spdf1 = filter( r-> !ismissing(r.Days) , spdf1 )
spdf1.p .= .1

spdf = DataFrame( Days = vec(spstats25')
	  	 , Count = repeat( [Symbol("1 case"), Symbol("3 cases"), :Infections], inner=size(spstats25,2)))
spdf2 = spdf[ spdf.Count .!= :Infections, :]
spdf2 = filter( r-> !ismissing(r.Days) , spdf2 )
spdf2.p .= 0.25

spdf = vcat( spdf1, spdf2 )

PlotlyJS.plot( cmplot(spdf;xcol=[:Count,:p],ycol=:Days,colorshift=1)... )

@show summarystats( spstats[3,:] )
#= Summary Stats:
Length:         1000
Missing Count:  96
Mean:           997.555310
Std. Deviation: 938.815997
Minimum:        1.000000
1st Quartile:   330.500000
Median:         751.000000
3rd Quartile:   1382.500000
Maximum:        7761.000000 =#


