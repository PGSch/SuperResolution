module SuperResolution
using Statistics
#cd("/Users/Patrick/Documents/GitHub/"
include("clusGAP.jl")
include("clusSR.jl")
include("evalSR.jl")
include("fctSR.jl")
include("importSR.jl")
include("plotSR.jl")
include("segmentSR.jl")
include("modelSR.jl")

greet()=print("Super Resolution - Blinking Artifacts")
filename1="locs_gamma_conf.csv"
filename2="labels_gamma.csv"
filepath1 = joinpath(@__DIR__, filename1)
filepath2 = joinpath(@__DIR__, filename2)
pStr1=filepath1
pStr2=filepath2
# pStr1="/Users/Patrick/Documents/GitHub/SuperResolution/src/locs_gamma_conf.csv"
# pStr2="/Users/Patrick/Documents/GitHub/SuperResolution/src/labels_gamma.csv"
global DataSim, RealLocs=importSR(pStr1,pStr2)
# using SuperResolution
SEG_SPLIT=100;
tmp1=SuperResolution.SegmentSR(SEG_SPLIT,SuperResolution.DataSim,SuperResolution.RealLocs); #segment data

iter_count=10;
SEG_RESTR=rand(tmp1.idx,iter_count);
EMP1=SuperResolution.EvalSR(tmp1,SEG_RESTR);
#
# EMP2=EMP1;
# #tmp2=tmp1;
# #
# asdf=map(i->sum(tmp1.N[CartesianIndex(tmp1.idx[i])]),1:iter_count)
# asdf2=map(i->(EMP1.ground_truth_K-EMP1.diff_WK_max_idx)[i]/EMP1.ground_truth_K[i],1:iter_count)
# plot([EMP1.ground_truth_K-EMP1.diff_WK_max_idx,log.(asdf),asdf2])
# plot(asdf2)
# #
##
# tmp2=SuperResolution.ModelSR(tmp1,(8,8),1)  #for segment create all models  WKLD
# tmp22=SuperResolution.ModelSR(tmp1,(8,8),2)  #for segment create all models KLD
# tmp3=SuperResolution.DistanceSR(tmp2,1)    #cluster WKLD
# tmp33=SuperResolution.DistanceSR(tmp2,2)   #cluster KLD
# tmp4=SuperResolution.clusGAP(tmp2,tmp22,1)  #Gap statistic
# tmp44=SuperResolution.clusGAP(tmp2,tmp22,2)  #Gap statistic
# #
# using SuperResolution
# tmp1=SuperResolution.SegmentSR(10,SuperResolution.DataSim,SuperResolution.RealLocs) #segment data
# tmp2=SuperResolution.ModelSR(tmp1,(2,3),1)  #for segment create all models  WKLD
# tmp22=SuperResolution.ModelSR(tmp1,(2,3),2)  #for segment create all models KLD
# tmp3=SuperResolution.DistanceSR(tmp2,1)    #cluster WKLD
# tmp33=SuperResolution.DistanceSR(tmp2,2)   #cluster KLD
# tmp4=SuperResolution.clusGAP(tmp2,tmp22,1)  #Gap statistic
# tmp44=SuperResolution.clusGAP(tmp2,tmp22,2)  #Gap statistic
# #
# id1=1
# id2=9
# tmp2=SuperResolution.ModelSR(tmp1,(id1,id2),1)  #for segment create all models  WKLD
# tmp22=SuperResolution.ModelSR(tmp1,(id1,id2),2)  #for segment create all models KLD
# tmp3=SuperResolution.DistanceSR(tmp2,1)    #cluster WKLD
# tmp33=SuperResolution.DistanceSR(tmp2,2)   #cluster KLD
# tmp4=SuperResolution.clusGAP(tmp2,tmp22,1)  #Gap statistic
# tmp44=SuperResolution.clusGAP(tmp2,tmp22,2)  #Gap statistic
#
# asdf=tmp4.W_K[1:(end-1)]-tmp4.W_K[2:end]
# plot(asdf)
#
# #
# scatter([tmp1.MU[2,3][:,1],tmp2.loc[end-5][:,1]],[tmp1.MU[2,3][:,2],tmp2.loc[end-5][:,2]])
# scatter(tmp1.MU[2,3][:,1],tmp1.MU[2,3][:,2])
# scatter(tmp2.loc[end-5][:,1],tmp2.loc[end-5][:,2])
# asdf=tmp4.W_K[1:(end-1)]-tmp4.W_K[2:end]
# plot(asdf)
#             K_SR::Int64
#             N_SR::Int64
#             val_bool::BitArray{1}
#             val_bool2::Array{BitArray{1},1}
#             singleDist_SR::Array{Array{Float64,2},1}
#             loc_SR::Array{Array{Array,1},1}
#             frame_SR::Array{Array{Int64,1},1}
#             cov_SR::Array{Array{Array{Float64,2},1},1}
#             D_SR::Array{Float64,1}
#             W_K::Float64
#             weight_SR::Array{Array{Float64,2},1}
# tmp1=SuperResolution.SegmentSR(10,SuperResolution.DataSim,SuperResolution.RealLocs)
# SAMPLE1=SuperResolution.ctrans(rand(tmp2.model[end-1],60))
# gather_X1=hcat(tmp1.loc[8,8][:,2],SAMPLE1[:,1])
# scatter(SuperResolution.𝕐sim[:,1],SuperResolution.𝕐sim[:,2])
# scatter(tmp1.loc[8,8][:,1],tmp1.loc[8,8][:,2])
# tmp2=SuperResolution.ModelSR(tmp1,(8,8))
# tmp3=(map(i->mean(tmp2.weight[i]),2:(length(tmp2.weight)-1)))
# tmp4=(map(i->mean(tmp2.param[i]),1:(length(tmp2.param)-1)))
# plot(tmp3)
# plot(tmp4)
# scatter(tmp1.loc[8,8][:,1],tmp1.loc[8,8][:,2])
# plot(contourf(exp(abs.(tmp2.wkld[end-8]-tmp3.wkld[end-8]))))
# global NN=20	# partition of space in (NN-1)^2 parts
#
# global evalErrSINGLE=[[0.0],]
# global evalErrGLOBAL=[[0.0],]
# global evalmeanErrGLOBAL=[0.0,]
# global evalcountErrGLOBAL=[0.0,]
# global evalcountErrGLOBAL2=[0,]
# #global evalSEED=[0,]
#
# global partN_valid=[0,]
# global partN2_valid=[0,]
#
# global estim_mol_numberGLOBAL=[0,]
# global sim_mol_numberGLOBAL=[0,]
# global obsGLOBAL=[0,]
#
# global 𝕐estAll=Array{Float64,2}(undef, 0,2)
# global 𝕐simAll=Matrix{Float64}(DataSim[[4,5],])
#
# global 𝕄simAll= Vector{Int64}(DataSim[1,])
# global 𝛍simAll= Matrix{Float64}(undef,length(unique(𝕄simAll)),2)
# #For unique molecule numbers compute means (needed for evaluation)
# for i=1:length(unique(𝕄simAll))
# 	global 𝛍simAll[i,:]=mean(𝕐simAll[𝕄simAll .== unique(𝕄simAll)[i],1:2]',dims=2)[:]
# end
end
