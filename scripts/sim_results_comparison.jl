# Comparing simulation results using different parameters

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

using Plots

### Load different simulation results
# Original preliminary results
sims1 = load("covidlike-1.0-sims.jld2", "sims")
# Plus age disaggregation of contact degree and assortativity
sims2 = load("covidlike-1.2-sims.jld2", "sims") # Note that the filename order is different from the sims name order
# Plus age disaggregation of contact degree and assortativity, and age disaggregation of commuting (16+ years only)
sims3 = load("covidlike-1.1-sims.jld2", "sims")


### Cumulative infections
times = range( 0.0, 60.0, length = 200 ) 	
# sims1 manipulations
it1 = reduce(hcat,  map( sims1 ) do s1
	map( t-> sum( s1.G.tinf .< t ), times )
end 
)
bdata1 = mapslices( i->quantile(i, [.5, .025, .975] ), it1, dims=2) 
# sims2 manipulations
it2 = reduce(hcat,  map( sims2 ) do s2 
	map( t-> sum( s2.G.tinf .< t ), times )
end 
)
bdata2 = mapslices( i->quantile(i, [.5, .025, .975] ), it2, dims=2) 
# sims3 manipulations
it3 = reduce(hcat,  map( sims3 ) do s3 
	map( t-> sum( s3.G.tinf .< t ), times )
end 
)
bdata3 = mapslices( i->quantile(i, [.5, .025, .975] ), it3, dims=2) 
# Create plot
fband = Figure(resolution=(600, 300)#3.5*72, 3.25*72)
                , fontsize=11);
# Add or remove yscale depending on required plot type
ax = Axis(fband[1, 1], xlabel="Days after first introduction", ylabel="Cumulative infections", title=""
	#, yscale=log10
	  , limits = (nothing,(minimum(bdata1[bdata1[:,2].>0,2])*0.9, maximum(bdata1[:,3])*1.1))
	);
# Add plot for sims1
band!(ax, times, bdata1[:,2], bdata1[:,3], color=(:blue, .15));
lines!(ax, times, bdata1[:,1], color=:black, label="Preliminary results");
# Add plot for sims2
band!(ax, times, bdata2[:,2], bdata2[:,3], color=(:pink, .15));
lines!(ax, times, bdata2[:,1], color=:red, label="Age disaggregation (contact degree & assortativity)");
# Add plot for sims3
band!(ax, times, bdata3[:,2], bdata3[:,3], color=(:lightgreen, .15));
lines!(ax, times, bdata3[:,1], color=:green, label="Age disaggregation (cd & assortativity & commuting)");
# Add legend
#axislegend(ax; position=:lt) 
Legend(fband[1,2],ax)
fband
# Select depending on plot type
save( "covidlike-1.0_to_1.2-fband.pdf", fband ) 
save( "covidlike-1.0_to_1.2-fband-log10.pdf", fband ) 


###############################################
### Distributions of time to detection (TD)

# time to detect
# Create function to generate data for plotting
function td_for_plot(; sim_data, sample_proportion)
    sp = sample_proportion
    spstats = map(sim_data)  do s 
		ss = sampleforest( s , sp )
		t0 = ss.firstsample
		t3 = (ss.n>2) ? ss.tsample[3] : missing
		i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
	    [ t0, t3, i0 ]
    end
    spstats = reduce( hcat, spstats )
    return(spstats)
end

spstats_sims1_sp10 = td_for_plot(sim_data=sims1, sample_proportion=0.1)
spstats_sims1_sp50 = td_for_plot(sim_data=sims1, sample_proportion=0.5)
spstats_sims2_sp10 = td_for_plot(sim_data=sims2, sample_proportion=0.1)
spstats_sims2_sp50 = td_for_plot(sim_data=sims2, sample_proportion=0.5)
spstats_sims3_sp10 = td_for_plot(sim_data=sims3, sample_proportion=0.1)
spstats_sims3_sp50 = td_for_plot(sim_data=sims3, sample_proportion=0.5)

@show summarystats( spstats_sims1_sp10[3,:] )
@show summarystats( spstats_sims1_sp50[3,:] )

# Manipulate data for plotting
function generate_spdf(;spstats,sample_proportion, sims_n)
    spdf = DataFrame( Days = vec(spstats')
	              	 , Count = repeat( [Symbol("1 case"), Symbol("3 cases"), :Infections], inner=size(spstats,2)))
    spdf1 = spdf[ spdf.Count .!= :Infections, :]
    spdf1 = filter( r-> !ismissing(r.Days) , spdf1 )
    spdf1.p .= sample_proportion
    spdf1.sims_n .= sims_n
    return( spdf1 )
end

spdf_sims1_sp10 = generate_spdf( spstats = spstats_sims1_sp10, sample_proportion = 0.1, sims_n = "sim1")
spdf_sims1_sp50 = generate_spdf( spstats = spstats_sims1_sp50, sample_proportion = 0.5, sims_n = "sim1")
spdf_sims2_sp10 = generate_spdf( spstats = spstats_sims2_sp10, sample_proportion = 0.1, sims_n = "sim2")
spdf_sims2_sp50 = generate_spdf( spstats = spstats_sims2_sp50, sample_proportion = 0.5, sims_n = "sim2")
spdf_sims3_sp10 = generate_spdf( spstats = spstats_sims3_sp10, sample_proportion = 0.1, sims_n = "sim3")
spdf_sims3_sp50 = generate_spdf( spstats = spstats_sims3_sp50, sample_proportion = 0.5, sims_n = "sim3")

#spdf = vcat( spdf1, spdf2 )
spdf_sims1 = vcat( spdf_sims1_sp10, spdf_sims1_sp50 )
PlotlyJS.plot( cmplot(spdf_sims1;xcol=[:Count,:p],ycol=:Days,colorshift=1)... )
spdf_sims2 = vcat( spdf_sims2_sp10, spdf_sims2_sp50 )
PlotlyJS.plot( cmplot(spdf_sims2;xcol=[:Count,:p,:sims_n],ycol=:Days,colorshift=2)... )
spdf_sims3 = vcat( spdf_sims3_sp10, spdf_sims3_sp50 )
PlotlyJS.plot( cmplot(spdf_sims3;xcol=[:Count,:p,:sims_n],ycol=:Days,colorshift=5)... )
# All together
spdf_sims_all = vcat( spdf_sims1_sp10, spdf_sims1_sp50
                     ,spdf_sims2_sp10, spdf_sims2_sp50
                     ,spdf_sims3_sp10, spdf_sims3_sp50 )
# Plot
PlotlyJS.plot( cmplot(spdf_sims_all;xcol=[:Count,:p,:sims_n],ycol=:Days,colorshift=1)... )
# Plot in order of number of cases
PlotlyJS.plot( cmplot(sort!(spdf_sims_all,[:Count,:p,:sims_n]);xcol=[:Count,:p,:sims_n],ycol=:Days,colorshift=1)... )

# With different colouring (not working)
unique_sims = unique(spdf_sims_all.sims_n)
sims_palette = palette(:viridis, length(unique_sims))  # or any other palette

# Map each sims_n to a color
sims_color_map = Dict(s => sims_palette[i] for (i, s) in enumerate(unique_sims))

# Then pass colors to CMPlot or PlotlyJS accordingly, e.g.,
sims_colors = [sims_color_map[s] for s in spdf_sims_all.sims_n]

# Plot and save file
fig_td_dist = PlotlyJS.plot( cmplot( sort!(spdf_sims_all,[:Count,:p,:sims_n])
                            ;xcol=[:Count,:p,:sims_n]
                            ,ycol=:Days
                            ,colorshift=1
                            #,color = sims_colors
                            )... )
fig_td_dist
save( "covidlike-1.0_to_1.2-td_dist.pdf", fig_td_dist ) 


### Extract data

# Extract mean time to detection (TD) in days for 1 and 3 cases (TD3)
# Function to return information
# TD at sampling proportion = 10%

function td_summary(; spdf, bdata, sample_proportion, times=range(0.0,60.0,length=200) )
    
    # TD
    spdf_1case = filter(row -> row.Count ==(Symbol("1 case")) && row.p == sample_proportion, spdf)
    td_mean = mean( spdf_1case.Days )
    td_lower = quantile( spdf_1case.Days,0.25 )
    td_upper = quantile( spdf_1case.Days,0.75 )
    # TD3
    spdf_3case = filter(row -> row.Count ==(Symbol("3 cases")) && row.p == sample_proportion, spdf)
    td3_mean = mean( spdf_3case.Days )
    td3_lower = quantile(spdf_3case.Days,0.25)
    td3_upper = quantile(spdf_3case.Days,0.75)
    # Cumulative infections
    # Indices for cumulative infections at time of detection (TD)
    idx_td_mean = minimum( findall(>(td_mean), times) )
    idx_td_lower = minimum( findall(>(td_lower), times) )
    idx_td_upper = minimum( findall(>(td_upper), times) )
    idx_td3_mean = minimum( findall(>(td3_mean), times) )
    idx_td3_lower = minimum( findall(>(td3_lower), times) )
    #idx_td3_upper = minimum( findall(>(td3_upper), times) )
    idx_td3_upper = isempty( findall(>(td3_upper), times)) ? length(times) : minimum( findall(>(td3_upper), times) )
    
    # Print summary
    # TD
    println("TD 1 case")
    println( "Mean TD with p=",sample_proportion," is ", round( td_mean, digits=2))
    println( "with IQR: ", round(td_lower,digits=2)," to ", round(td_upper,digits=2)," days")
    println("and the cumulative number of infections at time = ", round(times[idx_td_mean],digits=2)
            , " is ", bdata[idx_td_mean,1])
    println( "with IQR: ",bdata[idx_td_lower,1]," at ",round(times[idx_td_lower],digits=2)
            ,"and ",bdata[idx_td_upper,1]," at ",round(times[idx_td_upper],digits=2))
    println("Mean number of infections at TD (across all simulation) with p=",sample_proportion
            ," is ", round( mean( bdata[:,1] ), digits=2))
    
    # TD3
    println("TD 3 cases")
    println( "Mean TD3 with p=",sample_proportion," is ", round( td3_mean, digits=2))
    println( "with IQR: ", round(td3_lower,digits=2)," to ", round(td3_upper,digits=2)," days")
    println("and the cumulative number of infections at time = ", round(times[idx_td3_mean], digits=2)
            , " is ", bdata[idx_td3_mean,1])
    println( "with IQR: ",bdata[idx_td3_lower,1]," at ",round(times[idx_td3_lower],digits=2)," and ",bdata[idx_td3_upper,1]," at ",round(times[idx_td3_upper],digits=2))
    
    # Comparison
    println("Comparison between TD and TD3")
    println("TD3 is ", round(td3_mean / td_mean, digits=2)," x TD" )
    println("Cumulative Infections: TD3 is ", round(bdata[idx_td3_mean,1] / bdata[idx_td_mean,1], digits=2)," x TD" )

    
end

# sims1
td_summary( spdf = spdf_sims1_sp10, bdata = bdata1, sample_proportion = 0.10)
td_summary( spdf = spdf_sims1_sp50, bdata = bdata1, sample_proportion = 0.50)
# sims2
td_summary( spdf = spdf_sims2_sp10, bdata = bdata2, sample_proportion = 0.10)
td_summary( spdf = spdf_sims2_sp50, bdata = bdata2, sample_proportion = 0.50)
# sims3
td_summary( spdf = spdf_sims3_sp10, bdata = bdata3, sample_proportion = 0.10)
td_summary( spdf = spdf_sims3_sp50, bdata = bdata3, sample_proportion = 0.50)

# Mean number of cumulative infections at TD (1 case)
function mean_cumulative_infections(; sim=sims1)
    sps = (0.1,0.2,0.3,0.4,0.5) # sample proportions
    Y = map(sps) do sp 
	    y = map(sim)  do s 
		    ss = sampleforest( s , sp )
		    t0 = ss.firstsample
		    t3 = (ss.n>2) ? ss.tsample[3] : missing
		    i0 = (ss.n>0) ? sum( s.G.tinf .<= t0  ) : missing 
		    [ t0, t3, i0 ]
	    end
	    reduce( hcat, y )
    end
    Y10 = Y[1]; Y20 = Y[2]; Y30 = Y[3]; Y40 = Y[4]; Y50 = Y[5]
    Y10_filtered = filter(z->!ismissing(z), Y10[3,:] )
    Y20_filtered = filter(z->!ismissing(z), Y20[3,:] )
    Y30_filtered = filter(z->!ismissing(z), Y30[3,:] )
    Y40_filtered = filter(z->!ismissing(z), Y40[3,:] )
    Y50_filtered = filter(z->!ismissing(z), Y50[3,:] )
    println("Mean number of cumulative infections at TD:")# for $(string(sim)):")
    println("with p=",sps[1]," is ", mean(Y10_filtered))
    println("with p=",sps[2]," is ", mean(Y20_filtered),", a change of ", 1-mean(Y20_filtered)/mean(Y10_filtered))
    println("with p=",sps[3]," is ", mean(Y30_filtered),", a change of ", 1-mean(Y30_filtered)/mean(Y20_filtered))
    println("with p=",sps[4]," is ", mean(Y40_filtered),", a change of ", 1-mean(Y40_filtered)/mean(Y30_filtered))
    println("with p=",sps[5]," is ", mean(Y50_filtered),", a change of ", 1-mean(Y50_filtered)/mean(Y40_filtered))
end

mean_cumulative_infections( sim = sims1 )
mean_cumulative_infections( sim = sims2 )
mean_cumulative_infections( sim = sims3 )


###############################################
### Sensitivty to sample proportion
# t & i ~ p 
st0 = time() 
sps = range(.05, .5, length=20) |> collect 
# sims1
Y1 = map(sps) do sp1 
	y1 = map(sims1)  do s1 
		ss1 = sampleforest( s1 , sp1 )
		t0_1 = ss1.firstsample
		t3_1 = (ss1.n>2) ? ss1.tsample[3] : missing
		i0_1 = (ss1.n>0) ? sum( s1.G.tinf .<= t0_1  ) : missing 
		[ t0_1, t3_1, i0_1 ]
	end
	reduce( hcat, y1 )
end 
st1 = time() 
@show st1 - st0
Y1 = cat( Y1..., dims = 3 )
# sims2
Y2 = map(sps) do sp2 
	y2 = map(sims2)  do s2 
		ss2 = sampleforest( s2 , sp2 )
		t0_2 = ss2.firstsample
		t3_2 = (ss2.n>2) ? ss2.tsample[3] : missing
		i0_2 = (ss2.n>0) ? sum( s2.G.tinf .<= t0_2  ) : missing 
		[ t0_2, t3_2, i0_2 ]
	end
	reduce( hcat, y2 )
end 
st1 = time() 
@show st1 - st0
Y2 = cat( Y2..., dims = 3 )
# sims3
Y3 = map(sps) do sp3 
	y3 = map(sims3)  do s3 
		ss3 = sampleforest( s3 , sp3 )
		t0_3 = ss3.firstsample
		t3_3 = (ss3.n>2) ? ss3.tsample[3] : missing
		i0_3 = (ss3.n>0) ? sum( s3.G.tinf .<= t0_3  ) : missing 
		[ t0_3, t3_3, i0_3 ]
	end
	reduce( hcat, y3 )
end 
st1 = time() 
@show st1 - st0
Y3 = cat( Y3..., dims = 3 )

#f = Figure(size=(350,600));
f = Figure(size=(800,350));
rowtits = ["Day 1st detection", "Day 3rd detection", "Cumulative infections at detection"] 
for ri in 1:3
	#ri = 2
	#a = Axis( f[ri,1] , title = rowtits[ri], xlabel = "ICU sample proportion", ylabel="", yscale=(ri == 3 ? log10 : identity)) ;
    a = Axis( f[1,ri] , title = rowtits[ri], xlabel = "ICU sample proportion", ylabel="", yscale=(ri == 3 ? log10 : identity)) ;
	x = repeat(sps,inner=1000)
	y1 = vec(Y1[ri,:,:]) # sims1
    y2 = vec(Y2[ri,:,:]) # sims2
    y3 = vec(Y3[ri,:,:]) # sims3
	y_filtered = filter(z->!ismissing(z), y1 )
	#ylims!( a, quantile(y_filtered, .01), quantile(y_filtered,.99))
	Makie.ylims!( a, quantile(y_filtered, .01), quantile(y_filtered,.99))
	bw = diff( sps[1:2] )[1]/2 
	Makie.scatter!(a, x .+ bw.*(rand(length(x)).-.5), y1, color = :blue, alpha = .3, markersize=1 );
    Makie.scatter!(a, x .+ bw.*(rand(length(x)).-.5), y2, color = :red, alpha = .3, markersize=1 );
    Makie.scatter!(a, x .+ bw.*(rand(length(x)).-.5), y3, color = :green, alpha = .3, markersize=1 );
	if ri==3  
		m1 = lm(@formula(y ~ x), DataFrame(x = x, y = log.(y1)))
        m2 = lm(@formula(y ~ x), DataFrame(x = x, y = log.(y2)))
        m3 = lm(@formula(y ~ x), DataFrame(x = x, y = log.(y3)))
		py1 = exp.( predict( m1, DataFrame(x=sps)) )
        py2 = exp.( predict( m2, DataFrame(x=sps)) )
        py3 = exp.( predict( m3, DataFrame(x=sps)) )
		#exp.( predict( m, DataFrame(x=0.4)) )
		lines!(a, sps, py1 , color=:blue, linewidth=2);
        lines!(a, sps, py2 , color=:red, linewidth=2);
        lines!(a, sps, py3 , color=:green, linewidth=2);
        # Print model values for number of infections at sample proportion = {0.1,0.5}
        println("Cumulative infections at TD with sample proportion = 0.1 ")
        println("Cumulative infections at TD with sample proportion = 0.1 ")
	else 
		m1 = lm(@formula(y ~ x), DataFrame(x = x, y = y1))
        m2 = lm(@formula(y ~ x), DataFrame(x = x, y = y2))
        m3 = lm(@formula(y ~ x), DataFrame(x = x, y = y3))
		py1 = predict( m1, DataFrame(x=sps)) 
        py2 = predict( m2, DataFrame(x=sps)) 
        py3 = predict( m3, DataFrame(x=sps)) 
		lines!(a, sps, py1 , color=:blue, linewidth=2);
        lines!(a, sps, py2 , color=:red, linewidth=2);
        lines!(a, sps, py3 , color=:green, linewidth=2);
	end
	@show rowtits[ri] 
	@show DataFrame( sprop = sps, var = py1)#, var1 = py2, var1 = py3)
end 
f
save( "covidlike-1.0_to_1.2-sampleplots-ti_p.pdf", f )

