using SuperResolution


mutable struct EvalSR
    _selectIDX::Array{Tuple{Int64,Int64},1}
    wModel::Array{SuperResolution.ModelSR,1}
    RefModel::Array{SuperResolution.ModelSR,1}
    distSR::Array{SuperResolution.DistanceSR,1}
    gapSR::Array{SuperResolution.clusGAP,1}
    diff_WK_max_idx::Array{Int64,1}
    ground_truth_K::Array{Int64,1}
    k_corr_rate::Float64
    function EvalSR(SEG::SuperResolution.SegmentSR,RESTR_SEG::Array{Tuple{Int64,Int64},1})#RESTR_SEG::Int64)
        obj=new()
        obj._selectIDX=RESTR_SEG
        NUM_SEG=length(RESTR_SEG)
        #obj._selectIDX=SEG.idx[1:RESTR_SEG]
        obj.wModel=Array{SuperResolution.ModelSR,1}(undef,NUM_SEG)
        obj.RefModel=Array{SuperResolution.ModelSR,1}(undef,NUM_SEG)
        obj.distSR=Array{SuperResolution.DistanceSR,1}(undef,NUM_SEG)
        obj.gapSR=Array{SuperResolution.clusGAP,1}(undef,NUM_SEG)
        obj.diff_WK_max_idx=Array{Int64,1}(undef,NUM_SEG)
        obj.ground_truth_K=Array{Int64,1}(undef,NUM_SEG)

        @time(for i=1:NUM_SEG
                print("################################################\n")
                print("#..start segment $(SEG.idx[i])-$(i)#####################\n")
                print("################################################\n")
                obj.wModel[i]=SuperResolution.ModelSR(SEG,SEG.idx[i],1)  #for segment create all models  WKLD
                obj.RefModel[i]=SuperResolution.ModelSR(SEG,SEG.idx[i],2)  #for segment create all models KLD
                obj.distSR[i]=SuperResolution.DistanceSR(obj.wModel[i],1)    #cluster WKLD
                #tmp33=SuperResolution.DistanceSR(obj.wModel,2)   #cluster KLD
                obj.gapSR[i]=SuperResolution.clusGAP(obj.wModel[i],obj.RefModel[i],1)  #Gap statistic
                #tmp44=SuperResolution.clusGAP(obj.wModel,obj.RefModel,2)  #Gap statistic
        end)
        obj.diff_WK_max_idx=map(i->obj.gapSR[i].diff_WK_max_idx[1],1:length(obj.gapSR))
        obj.ground_truth_K=map(i->SEG.K[CartesianIndex(SEG.idx[i])],1:NUM_SEG)
        obj.k_corr_rate=sum(abs.(obj.ground_truth_K-obj.diff_WK_max_idx))/sum(obj.ground_truth_K)
        print("################################################\n")
        print("#correct detection rate $(obj.k_corr_rate)#####################\n")
        print("################################################\n")
        return obj
    end
end
export EvalSR
