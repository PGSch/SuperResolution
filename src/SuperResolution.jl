module SuperResolution
using Statistics
#cd("/Users/Patrick/Documents/GitHub/"
# include("clusGAP.jl")
# include("clusSR.jl")
# include("evalSR.jl")
# include("fctSR.jl")
# include("importSR.jl")
# include("plotSR.jl")
# include("segmentSR.jl")
# include("wkldSR.jl")

greet()=print("Super Resolution - Blinking Artefacts")
# Write your package code here.

pStr1="/Users/Patrick/Documents/GitHub/SuperResolution/src/locs_gamma_conf.csv"
pStr2="/Users/Patrick/Documents/GitHub/SuperResolution/src/labels_gamma.csv"
global DataSim, RealLocs=importSR(pStr1,pStr2)

global NN=20	# partition of space in (NN-1)^2 parts

global evalErrSINGLE=[[0.0],]
global evalErrGLOBAL=[[0.0],]
global evalmeanErrGLOBAL=[0.0,]
global evalcountErrGLOBAL=[0.0,]
global evalcountErrGLOBAL2=[0,]
#global evalSEED=[0,]

global partN_valid=[0,]
global partN2_valid=[0,]

global estim_mol_numberGLOBAL=[0,]
global sim_mol_numberGLOBAL=[0,]
global obsGLOBAL=[0,]

global 𝕐estAll=Array{Float64,2}(undef, 0,2)
global 𝕐simAll=Matrix{Float64}(DataSim[[4,5],])

global 𝕄simAll= Vector{Int64}(DataSim[1,])
global 𝛍simAll= Matrix{Float64}(undef,length(unique(𝕄simAll)),2)
#For unique molecule numbers compute means (needed for evaluation)
for i=1:length(unique(𝕄simAll))
	global 𝛍simAll[i,:]=mean(𝕐simAll[𝕄simAll .== unique(𝕄simAll)[i],1:2]',dims=2)[:]
end
end
