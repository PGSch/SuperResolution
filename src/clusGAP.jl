using SuperResolution
include("clusSR.jl")

mutable struct clusGAP
    modus::Int64
    w_KLD_Distance_SR::SuperResolution.DistanceSR
    KLD_Distance_SR::SuperResolution.DistanceSR
    W_K::Array{Float64,1}
    W_K_star::Array{Float64,1}
    Gap_n::Array{Float64,1}
    Gap::Array{Float64,1}
    s_k::Array{Float64,1}
    sd_k::Array{Float64,1}
    k_hat::Int64
    function clusGAP(ModelWKLD_SR::SuperResolution.ModelSR,ModelRef_SR::SuperResolution.ModelSR,modINPUT=1)
        obj=new()
        obj.modus=modINPUT
        obj.w_KLD_Distance_SR=SuperResolution.DistanceSR(ModelWKLD_SR,1)
        obj.W_K=map(i->obj.w_KLD_Distance_SR.Distance[i].W_K,1:length(obj.w_KLD_Distance_SR.Distance))
        if obj.modus==1
            obj.KLD_Distance_SR=SuperResolution.DistanceSR(ModelRef_SR,2)
            obj.W_K_star=map(i->obj.KLD_Distance_SR.Distance[i].W_K,1:length(obj.KLD_Distance_SR.Distance))
            obj.Gap_n=map(i->log(obj.KLD_Distance_SR.Distance[i].W_K)-log(obj.w_KLD_Distance_SR.Distance[i].W_K),1:(length(obj.w_KLD_Distance_SR.Distance)-1))
        elseif obj.modus==2
            obj.KLD_Distance_SR=SuperResolution.DistanceSR(ModelWKLD_SR,2)
            obj.W_K_star=map(i->obj.KLD_Distance_SR.Distance[i].W_K,1:length(obj.KLD_Distance_SR.Distance))
            obj.Gap_n=map(i->log(obj.KLD_Distance_SR.Distance[i].W_K)-log(obj.w_KLD_Distance_SR.Distance[i].W_K),1:(length(obj.w_KLD_Distance_SR.Distance)-1))
        end
        obj
    end
end
export clusGAP
