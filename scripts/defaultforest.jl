using NBMscape 
using CSV 

fo = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100)
CSV.write("defaultforest-D.csv", fo.D )
CSV.write("defaultforest-G.csv", fo.G )


fo1 = simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=4)
D = fo1.D[ .!isnothing.( fo1.D[:,1]),: ]
CSV.write("defaultforest-gen4-D.csv",D )
CSV.write("defaultforest-gen4-G.csv", fo1.G )

