using NBPMscape 
using CSV 
using Revise

fo = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100)
CSV.write("defaultforest-D.csv", fo.D )
CSV.write("defaultforest-G.csv", fo.G )


fo1 = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=4)
D = fo1.D[ .!isnothing.( fo1.D[:,1]),: ]
CSV.write("defaultforest-gen4-D.csv",D )
CSV.write("defaultforest-gen4-G.csv", fo1.G )

P =  merge( NBPMscape.P , (;infectivity=3)); 
fo2 = simforest(P; initialtime=0.0, maxtime=60.0, maxgenerations=4)
D = fo2.D[ .!isnothing.( fo2.D[:,1]),: ]
CSV.write("defaultforest-fo2-gen4-D.csv",D )
CSV.write("defaultforest-fo2-gen4-G.csv", fo2.G )


P =  merge( NBPMscape.P , (;infectivity=3)); 
fo3 = simforest(P; initialtime=0.0, maxtime=60.0, maxgenerations=6)
D = fo3.D[ .!isnothing.( fo3.D[:,1]),: ]
CSV.write("defaultforest-fo3-D.csv",D )
CSV.write("defaultforest-fo3-G.csv", fo3.G )


