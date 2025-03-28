using Revise
using NBPMscape 
using JLD2 
using StatsBase
using CairoMakie
using CMPlot 
using PlotlyJS
using DataFrames
using GLM, Statistics
using KernelDensity


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
sp = .1 
spstats = map(sims)  do s 
		ss = sampleforest( s , sp )
		t0 = ss.firstsample
		t3 = (ss.n>2) ? ss.tsample[3] : missing
		i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
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

#= Length:         1000
Missing Count:  94
Mean:           1012.248344
Std. Deviation: 984.907957
Minimum:        6.000000
1st Quartile:   321.750000
Median:         723.000000
3rd Quartile:   1403.750000
Maximum:        9760.000000
=#

# t & i ~ p 

st0 = time() 
sps = range(.05, .5, length=20) |> collect 
Y = map(sps) do sp 
	y = map(sims)  do s 
		ss = sampleforest( s , sp )
		t0 = ss.firstsample
		t3 = (ss.n>2) ? ss.tsample[3] : missing
		i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
		[ t0, t3, i0 ]
	end
	reduce( hcat, y )
end 
st1 = time() 
@show st1 - st0
Y = cat( Y..., dims = 3 )

f = Figure(size=(350,600));
rowtits = ["Day 1st detection", "Day 3rd detection", "Cumulative infections at detection"] 
for ri in 1:3 
	a = Axis( f[ri,1] , title = rowtits[ri], xlabel = "ICU sample proportion", ylabel="", yscale=(ri == 3 ? log10 : identity)) ;
	x = repeat(sps,inner=1000)
	y = vec(Y[ri,:,:]) 
	y1 = filter(z->!ismissing(z), y )
	ylims!( a, quantile(y1, .01), quantile(y1,.99))
	bw = diff( sps[1:2] )[1]/2
	scatter!(a, x .+ bw.*(rand(length(x)).-.5), y, color = :blue, alpha = .3, markersize=1 );
	if ri==3  
		m = lm(@formula(y ~ x), DataFrame(x = x, y = log.(y)))
		py = exp.( predict( m, DataFrame(x=sps)) )
		lines!(a, sps, py , color=:red, linewidth=2);
	else 
		m = lm(@formula(y ~ x), DataFrame(x = x, y = y))
		py = predict( m, DataFrame(x=sps)) 
		lines!(a, sps, py , color=:red, linewidth=2);
	end
	@show rowtits[ri] 
	@show DataFrame( sprop = sps, var = py)
end 
f
save( "covidlike-1.0-sampleplots-ti_p.pdf", f )
