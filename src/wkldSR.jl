using SuperResolution
using Distributions
include("fctSR.jl")

###4.. Setup KLClustering

##	Generate all possible clustering models for a _specific segment_
#	For final use, create such ModelSR for every valid SEG.idx
mutable struct ModelSR
		#Segment-restricted features
	dim::Int64
	loc::Array{Array{Float64,2},1}	#localisations of cluster populations
	model::Array{MixtureModel,1}	#MixtureModel
	member::Array{Array{Int64,1},1}	#cluster ID for all localisations
	clus_ID::Array{Array{Int64,1},1}	#cluster ID for remaining populations (needed for dim2WKLDfull; might replace with unique() and findall() for indices)

	cov::Array{Array{Array{Float64,2},1},1}	#covariance matrices
	comp_mix::Array{Array{Float64,1},1}	#convex coefficient of all MixtureModel components
	wkld::Array{Array{Float64,2},1}	#distance matrices KLD & comp_mix
	weight::Array{Array{Float64,2},1}	#distance matrices KLD & comp_mix
	geom::Array{Array{Geometric{Float64},1},1}
	param::Array{Array{Float64,1},1}

	K::Array{Int64,1}	#number of molecules in segment (according to ground_truth)
##Function to initialize structure
	function ModelSR(SEG::SuperResolution.SegmentSR,idx::Tuple{Int64,Int64}) #idx from list SEG.idx[j]
		obj=new()
		obj.dim=size(SEG.all_loc[:,1:(end-1)])[2]	#all_loc last col is framenumbers
		obj.loc=Array{Array{Float64,2},1}(undef,SEG.N[CartesianIndex(idx)]) #
		obj.member=Array{Array{Int64,1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.clus_ID=Array{Array{Int64,1},1}(undef,SEG.N[CartesianIndex(idx)])

		obj.model=Array{MixtureModel,1}(undef,SEG.N[CartesianIndex(idx)])
		obj.cov=Array{Array{Array{Float64,2},1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.comp_mix=Array{Array{Float64,1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.wkld=Array{Array{Float64,2},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.weight=Array{Array{Float64,2},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.geom=Array{Array{Geometric{Float64},1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.param=Array{Array{Float64,1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.K=Array{Int64,1}(undef,SEG.N[CartesianIndex(idx)])

		obj.member[1]=collect(1:SEG.N[CartesianIndex(idx)])#Array{Array{Int64,1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.clus_ID[1]=collect(1:SEG.N[CartesianIndex(idx)])#Array{Array{Int64,1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.cov[1]=repeat([1.0*eye(obj.dim)],SEG.N[CartesianIndex(idx)])#Array{Array{Array{Float64,2},1},1}(undef,SEG.N[CartesianIndex(idx)])
		obj.comp_mix[1]=1/SEG.N[CartesianIndex(idx)]*ones(1,SEG.N[CartesianIndex(idx)])[1,:]#Array{Array{Float64,2},1}(undef,SEG.N[CartesianIndex(idx)])

		sKLD=[dim2adj(SEG.loc[CartesianIndex(idx)]),]	#dim2adj: 452×3 Array{Float64,2} --> 452-element Array{Array,1}
		N=size(sKLD[1])[1]

		global i_inc=1
		sKLD[1],obj.comp_mix[1],obj.cov[1],obj.clus_ID[1],obj.member[1],obj.wkld[1],obj.weight[i_inc+1],obj.geom[1]=dim2WKLDfull(sKLD[i_inc],obj.cov[i_inc],obj.comp_mix[i_inc],SEG.frame[CartesianIndex(idx)],obj.clus_ID[i_inc],obj.member[i_inc])
		# sKLD[1],ωall[1],Σall[1],Aall[1],Acoll[1],KLDeval[1]=dim2WKLDfull(sKLD[i_inc],Σall[i_inc],ωall[i_inc],Aall[i_inc],Acoll[i_inc])
##1..Extract needed features
		global 𝕐simclAll=[SEG.loc[CartesianIndex(idx)],] #x,y coordinates of all observations [MANDATORY FOR MODEL SELECTION/EVALUATION]
		obj.model[1]=dim2mix(SEG.loc[CartesianIndex(idx)][:,1:2])
		obj.wkld[1]=deepcopy(SEG.loc[CartesianIndex(idx)])

		###4.. KLClustering
		print("number of observations = $(SEG.N[CartesianIndex(idx)])\n")
		@time(
		while size(obj.clus_ID[i_inc],1)>1
			print("i_inc= $(i_inc)\n")
			tmpsKLD,obj.comp_mix[i_inc+1],obj.cov[i_inc+1],obj.clus_ID[i_inc+1],obj.member[i_inc+1],obj.wkld[i_inc+1],obj.weight[i_inc+1],obj.geom[i_inc+1]=dim2WKLDfull(sKLD[i_inc],obj.cov[i_inc],obj.comp_mix[i_inc],SEG.frame[CartesianIndex(idx)],obj.clus_ID[i_inc],obj.member[i_inc])
			# tmpA=tmpA[map(i->tmpA[i]!=0,1:length(tmpA))]
			push!(sKLD,tmpsKLD)#x,y values + pdf value of Framenumber
			push!(𝕐simclAll,dim2adj(sKLD[i_inc+1]))#x,y values from sKLD variable
			push!(obj.model,dim2mix(dim2adj(sKLD[i_inc+1]),obj.comp_mix[i_inc+1]))#MixtureModel (PDF)
			global i_inc+=1
		end
		)
		for j=1:(length(obj.geom)-1)
			tete=map(i->params(obj.geom[j][i])[1],1:length(obj.geom[j]))
			obj.param[j]=tete
		end
		# for j=1:(length(tmp2.geom)-1)
		# 	tete=map(i->params(tmp2.geom[j][i])[1],1:length(tmp2.geom[j]))
		# 	asdf[j]=tete
		# end
		obj
	end
end

export ModelSR
##
# #x,y coordinates of cluster centers, for each iteration
# 𝕐simclAll=map(i->𝕐simclAll[i][:,1:2],1:length(𝕐simclAll))
#
# ###5.. Cluster Evaluation
# #sMMsim=dim2mix(𝛍sim,𝛚sim)#GMM-Distribution based on mean ("RealLoc") mol_IDs, for comparison
#
# print("##########################\n")
# print("Compute Model Selection...\n")
# print("##########################\n")
#
# xRangeSel=collect(1:(i_inc))
# Knum=map(i->size(𝕐simclAll[i],1),collect(1:size(𝕐simclAll,1))) #number of clusters for each interation
#
# #compute number, index and locations of unique molecules
# Acoll_unique=map(i->length(unique(Acoll[i])),1:length(Acoll))
# #print("Acoll_unique:\n"); show(Acoll_unique); print("\n")
# #compute min|max|mean KLD values
# #2:length(obj.wkld) because the 1st entry is [Y,] !!!!
# minKLDeval=map(i->minimum(obj.wkld[i][1e-10 .<obj.wkld[i]]),1:length(obj.wkld)-1)
# push!(minKLDeval,sort(obj.wkld[end-2][minKLDeval[end-1].<obj.wkld[end-2]])[1]) #pseudo minKLD value for 1 Cluster
#
# #compute absolute minKLD changes
# abs_minKLDeval=abs.(minKLDeval[1:end-1]-minKLDeval[2:end])
# push!(abs_minKLDeval,0)
# abs_minKLDeval[2:end]=abs_minKLDeval[1:end-1]
# abs_minKLDeval[1]=0
# #print("abs_minKLDeval:\n"); show(abs_minKLDeval); print("\n")
#
# #𝕐simclAll[end-solvIDXend] is the clustering solution
# global solvIDX=findall(abs_minKLDeval.==maximum(abs_minKLDeval))[1] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# global solvNumClust=size(𝕐simclAll[solvIDX],1)
# print("estimated cluster for partition $(partN)_$(partN2): $(solvNumClust)\n")
# # print("length(minKLDeval)=$(length(minKLDeval))\n")
# # print("length(abs_delta_minKLDeval)=$(length(abs_minKLDeval))\n")
