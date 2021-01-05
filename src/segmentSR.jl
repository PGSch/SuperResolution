using PDMats
include("fctSR.jl")

mutable struct BoundsSegment
	x::Array{Float64,1}
	y::Array{Float64,1}
	function BoundsSegment(x,y)
		obj=new()
		obj.x=x
		obj.y=y
		obj
	end
end

mutable struct SegmentSR
	select::Function
	idx::Array{Array{Int}}	#list of valid/verified segment indices
		#Segment-restricted features
	ground_truth::Array{Array{Float64}}
	loc::Array{Array{Float64}}
	frame::Array{Array{Int64}}
	mol_ID::Array{Array{Int64}}
	MU::Array{Array{Float64}}
	K::Array{Array{Int}}	#number of molecules in segment (according to ground_truth)
	N::Array{Array{Int}}	#number of observations in segment
	bounds::Array{BoundsSegment}	#typeof(BoundsSegment)=BoundsSegment

		#Overall (unrestricted) features
	all_ground_truth::Array{Float64}
	all_loc::Array{Float64}
	all_frame::Array{Int64}
	all_mol_ID::Array{Int64}
##Function to initialize structure
	function SegmentSR(NN::Int64,df1::DataFrame,df2::DataFrame)
		DataSim=df1,RealLoc=df2
		obj=new()
		obj.idx=Array{Array{Int}}(undef,0)
##Segment-restricted features
		obj.ground_truth=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	ℝsim
		obj.loc=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	𝕐sim
		obj.frame=Array{Array{Int64,1},2}(undef,NN-1,NN-1)	#	𝔽sim
		obj.mol_ID=Array{Array{Int64,1},2}(undef,NN-1,NN-1)	#	𝕄sim
		obj.MU=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	mean of members for each mol_ID
		obj.K=Array{Array{Int64,1},2}(undef,NN-1,NN-1)
		obj.N=Array{Array{Int64,1},2}(undef,NN-1,NN-1)
		obj.bounds=Array{BoundsSegment,2}(undef,NN-1,NN-1)
##Overall (unrestricted) features
		obj.all_ground_truth=Matrix{Float64}(RealLoc[[2,3,1],])	#Real Localisations
		obj.all_mol_ID=Matrix{Float64}(RealLoc[[2,3,1],])	#Real Localisations
		obj.all_frame=Vector{Int64}(DataSim[3,])	#Framenumber
		obj.all_mol_ID=Vector{Int64}(DataSim[1,])	#mol_ID
##1..Extract needed features
			#Mandatory column order 𝕐sim:
			#[:,1]~x[nm]|[:,2]~y[nm]|[:,3]~franumber
		global 𝕐sim = Matrix{Float64}(DataSim[[4,5,3],])		#Observations
			#Mandatory column order 𝕄sim:
			#[:,1]~x[nm]|[:,2]~y[nm]|[:,3]~moleculeID
		global ℝsim = Matrix{Float64}(RealLoc[[2,3,1],])	#Real Localisations
		global 𝔽sim = Vector{Int64}(DataSim[3,])	#Framenumber
		global 𝕄sim = Vector{Int64}(DataSim[1,])	#mol_ID
		global amountM=Vector{Int64}()
##2.. Boundries/Ticks
		xmin=0.0 #dim2range(𝕐sim)[1][1]
		ymin=0.0 #dim2range(𝕐sim)[2][1]
		xmax=ceil(maximum([dim2range(𝕐sim)[1][end],dim2range(𝕐sim)[2][end]]))
		ymax=xmax
			#Alternatively, set data partition based on evaluated maximum measurements and random "u"
		Xpart=Vector{Float64}()
		for k=0:(NN-1) push!(Xpart,(xmax-xmin)*k/NN) end #coordinate ticks of segmentation
		Ypart=Xpart
	# xmin=Xpart[partN]
	# xmax=Xpart[partN+1]
	# ymin=Ypart[partN2]
	# ymax=Ypart[partN2+1]
		for j=1:(NN-1) for k=1:(NN-1) obj.bounds[j,k].x=[Xpart[j],Xpart[j+1]], obj.bounds[j,k].y=[Ypart[k],Ypart[k+1]] end end

##3.. Restrict Area
		for partN=1:(NN-1)
			for partN2=1:(NN-1)
				print("VVVVVVVVVVVVVVVVVVV\n")
				print("partN = $(partN)\n")
				print("partN2 = $(partN2)\n")
				xmin=Xpart[partN]
				xmax=Xpart[partN+1]
				ymin=Ypart[partN2]
				ymax=Ypart[partN2+1]
				#Crop sim-data based on partition
				obj.loc[partN,partN2],val=dim2crop(𝕐sim,xmin,xmax,ymin,ymax)	#val is the boolean for indices in contraint area
				push!(amountM,size(𝕐sim)[1])#######################substitute for object field obj.N???#########################################################################################################################
				obj.N[partN,partN2]=size(𝕐sim)[1]
##4..Verify if segments are valid under constraints
				if size(𝕐sim)[1] > 5 && size(unique(𝕄sim[val]),1)>1
					push!(partN_valid,partN)
					push!(partN2_valid,partN2)
					#Crop Framenumber, mol_ID and Real Localisations based on cropped sim-data
					𝔽sim = 𝔽sim[val]
					obj.frame[partN,partN2]= 𝔽sim[val]
					𝕄sim= 𝕄sim[val]	#mol_ID
					obj.mol_ID[partN,partN2]= 𝕄sim[val]
					ℝsim,val2 = dim2crop(ℝsim,xmin,xmax,ymin,ymax)
					obj.ground_truth[partN,partN2]= dim2crop(ℝsim,xmin,xmax,ymin,ymax)

					#For validation purpose, determine amount of simulated molecules
					𝛍sim= Matrix{Float64}(undef,length(unique(𝕄sim)),2)
					#For unique molecule numbers compute means
					for i=1:length(unique(𝕄sim))
						𝛍sim[i,:]=mean(𝕐sim[𝕄sim .== unique(𝕄sim)[i],1:2]',dims=2)[:]
					end
					obj.MU[partN,partN2]=𝛍sim
				end
			end
		end
		obj
	end
end
export segmentSR
##


#
# ##
#
# 				###4.. Setup KLClustering
# 				#𝕐sim[:,3]=map(i->pdf(G,𝕐sim[i,3]),1:length(𝕐sim[:,3]))
# 				Y=dim2adj(𝕐sim)
# 				𝛆=minimum(Y)
# 				𝛅 = Int64(3)
# 				global 𝕐simclAll=[𝕐sim,] #x,y coordinates of all observations
# 				sKLD=[Y,]
# 				#KLDeval=Array{Array{Float64,2},1}(undef,length(Y))
# 				global KLDeval=deepcopy([𝕐sim,])
# 				N=size(sKLD[1])[1]
# 				global sMM = [dim2mix(𝕐sim[:,1:2]),]
# 				global Wall = [1/N*ones(1,N)[1,:],]
# 				global ω = [1/N*ones(1,N)[1,:],]
# 				global Σall = [repeat([1.0*eye(𝛅)],N),]
# 				global Aall = [collect(1:N),] #cluster indices: not merged indices
# 				global Acoll = [collect(1:N),] #cluster indices: indices each observation to track clustering progress; entries are quivalent to 𝕄sim
# 				global κtmp=size(unique(𝕄sim),1)
# 				global ι=1
#
# 				#KLDeval[1]=dim2WKLDfull(sKLD[ι],Σall[ι],Wall[ι],α,Aall[ι],Acoll[ι])[6]
# 				sKLD[1],Wall[1],Σall[1],Aall[1],Acoll[1],KLDeval[1]=dim2WKLDfull(sKLD[ι],Σall[ι],Wall[ι],Aall[ι],Acoll[ι])
# 				#print("Aall[$(ι)]=$(Aall[ι])\n")
# 				###4.. KLClustering
# 				print("number of observations = $(N)\n");push!(obsGLOBAL,N)
# 				#print("amountM = $(N)\n")
# 				#print("κtmp before while = $(κtmp)\n")
# 				@time(
# 				while size(Aall[ι],1)>0
# 					#print("o=",o,"\n")
# 					# N=maximum(size(sKLD[i]))
# 					#print("Aall[$(ι)]=$(Aall[ι])\n")
# 					tmpsKLD,tmpW,tmpS,tmpA,tmpA_coll,tmpKLDeval=dim2WKLDfull(sKLD[ι],Σall[ι],Wall[ι],Aall[ι],Acoll[ι])
# 					tmpA=tmpA[map(i->tmpA[i]!=0,1:length(tmpA))]
#
# 					push!(Wall,tmpW)#Weight Coefficient
# 					push!(Σall,tmpS)#COV Matrices
# 					push!(Aall,tmpA)#Cluster Indices
# 					push!(Acoll,tmpA_coll)#Cluster Indices without removing collapsed observations
# 					tmpsKLD=tmpsKLD[map(i->!isnan(tmpsKLD[i][1]),1:length(tmpsKLD))]
# 					tmpW=tmpW[map(i->!isnan(tmpW[i]),1:length(tmpW))]
# 					push!(sKLD,tmpsKLD)#x,y values + pdf value of Framenumber
# 					push!(KLDeval,tmpKLDeval)
# 					push!(ω,tmpW)#redundant!!! remove/merge variable!!
# 					push!(𝕐simclAll,dim2adj(sKLD[ι]))#x,y values from sKLD variable
# 					𝕐simcl=dim2adj(sKLD[ι])
# 					push!(sMM,dim2mix(𝕐simcl,ω[ι]))#MixtureModel (PDF)
# 					global κtmp=size(sKLD[ι],1)
# 					global ι+=1
# 				end
# 				)
# 				#x,y coordinates of cluster centers, for each iteration
# 				𝕐simclAll3D=𝕐simclAll
# 				𝕐simclAll=map(i->𝕐simclAll[i][:,1:2],1:length(𝕐simclAll))
#
# 			###5.. Cluster Evaluation
# 				#sMMsim=dim2mix(𝛍sim,𝛚sim)#GMM-Distribution based on mean ("RealLoc") mol_IDs, for comparison
#
# 				print("##########################\n")
# 				print("Compute wKLD Solution...\n")
# 				print("##########################\n")
