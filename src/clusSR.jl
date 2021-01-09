using SuperResolution
using Distributions
include("fctSR.jl")
include("modelSR.jl")

mutable struct SimSR
    sample::Array{Array{Float64,2},1} #each array entry corresponds to a sample
    x::Array{Float64,2} #gathered x values of all simulations (for plot purpose)
    y::Array{Float64,2} #gathered y values
    function SimSR(model::MixtureModel,simN,simIT)  #simN..number of samples drawn // #simIT..repetitions of sampling
         obj=new()
         obj.sample=Array{Array{Float64,2},1}(undef,simIT)
         obj.x=Array{Float64,2}(undef,simN,simIT)
         obj.y=Array{Float64,2}(undef,simN,simIT)
         for i=1:simIT obj.sample[i]=ctrans(rand(model,simN)) end
         obj.x[:,1]=obj.sample[1][:,1]
         obj.y[:,1]=obj.sample[1][:,2]
         map(i->obj.x[:,1:i]=hcat(obj.x[:,1:(i-1)],obj.sample[i][:,1]),2:length(obj.sample))
         map(i->obj.y[:,1:i]=hcat(obj.y[:,1:(i-1)],obj.sample[i][:,2]),2:length(obj.sample))
         obj
    end
end

#   !!for a selected model K_SR!! (i.e., there are assumably K_SR clusters in this model),
#   compute the inner distances and filter members for each clus_ID
mutable struct SegDistanceSR
    K_SR::Int64
    N_SR::Int64
    val_bool::BitArray{1}
    val_bool2::Array{BitArray{1},1}
    singleDist_SR::Array{Array{Float64,2},1}
    loc_SR::Array{Array{Array,1},1}
    frame_SR::Array{Array{Int64,1},1}
    cov_SR::Array{Array{Array{Float64,2},1},1}
    D_SR::Array{Float64,1}
    W_K::Float64
    weight_SR::Array{Array{Float64,2},1}
    function SegDistanceSR(modelSR::SuperResolution.ModelSR,K_SR::Int64,modINPUT::Int64)
        obj=new()
        obj.K_SR=K_SR
        obj.val_bool=length.(unique.(modelSR.member)).==K_SR
        obj.N_SR=length(modelSR.member[obj.val_bool])
        obj.D_SR=zeros(K_SR)

        obj.val_bool2=map(i->modelSR.member[obj.val_bool][1].==modelSR.clus_ID[obj.val_bool][1][i],1:length(modelSR.clus_ID[obj.val_bool][1])) #Array{BitArray{1},1}
        obj.loc_SR=map(i->dim2adj(modelSR.loc[1][obj.val_bool2[i],:]),1:length(obj.val_bool2))   #dim2adj(Array{Float64,2} of x,y loc) ->  Array{Array,1}
        obj.frame_SR=map(i->modelSR.frame[obj.val_bool2[i]],1:length(obj.val_bool2))
        obj.cov_SR=map(i->modelSR.cov[1][obj.val_bool2[i]],1:length(obj.val_bool2))
        #obj.weight_SR=ones(K_SR)#map(i->modelSR.weight[1][obj.val_bool2[i],:],1:length(obj.val_bool2))
        # dim2KLDfull(MU::Array{Array,1},S::Array{Array{Float64,2},1},modusKLD=2,weightParaKLD=NaN)
        if modINPUT==1
            obj.singleDist_SR=map(i->dim2KLDfull_nonClust(obj.loc_SR[i],obj.cov_SR[i],obj.frame_SR[i])[1],1:K_SR)
            obj.weight_SR=map(i->dim2KLDfull_nonClust(obj.loc_SR[i],obj.cov_SR[i],obj.frame_SR[i])[2],1:K_SR)
            obj.D_SR=map(i->sum(obj.singleDist_SR[i]),1:K_SR)
            obj.W_K=sum(i->sum(obj.D_SR[i]/(2*sum(obj.val_bool2[i]))),1:obj.K_SR)
            obj
        elseif modINPUT==2
        #dim2KLDfull_nonClust(MU::Array{Array,1},S::Array{Array{Float64,2},1},SEG_FrNm)
            obj.singleDist_SR=map(i->dim2KLDfull(obj.loc_SR[i],obj.cov_SR[i]),1:K_SR)
            obj.D_SR=map(i->sum(obj.singleDist_SR[i]),1:K_SR)
            obj.W_K=sum(i->sum(obj.D_SR[i]/(2*sum(obj.val_bool2[i]))),1:obj.K_SR)
            obj
        end
    end
end
mutable struct DistanceSR
    Distance::Array{SuperResolution.SegDistanceSR,1}
    function DistanceSR(model_SR::SuperResolution.ModelSR,modINPUT::Int64)
        obj=new()
        obj.Distance=Array{SuperResolution.SegDistanceSR,1}(undef,model_SR.K[1])
        if modINPUT==1
            print("##################################################################\n")
            print("Computing full DistanceSR struct for var.....WEIGHTED KLD!!!\n")
            @time(
                for j=1:model_SR.K[1]
                    obj.Distance[j]=SegDistanceSR(model_SR,j,modINPUT)
                end
            )
        elseif modINPUT==2
            print("##################################################################\n")
            print("Computing full DistanceSR struct for var.....NON-WEIGHTED KLD!!!\n")
            @time(
                for j=1:model_SR.K[1]
                    obj.Distance[j]=SegDistanceSR(model_SR,j,modINPUT)
                end
            )
        end
        obj
    end
end
export SimSR
export SegDistanceSR
export DistanceSR


#unique(tmp5.Distance[end-1].singleDist_SR)
