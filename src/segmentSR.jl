using PDMats
include("fctSR.jl")
#tmp1=SuperResolution.SegmentSR(100,SuperResolution.DataSim,SuperResolution.RealLocs)
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
	#EXAMPLE:
	#tmp1.ground_truth[3] returns the ground truth localisation of populations
	#tmp1.loc[3] returns all measurements [x,y,FrNm] for the 3rd segment
	#tmp1.MU[3] returns the means for each population in the 3rd segment
	select::Function
	idx::Array{Tuple{Int64,Int64},1}	#list of valid/verified segment indices
		#Segment-restricted features
	ground_truth::Array{Array{Float64,2},2}	#loc of populations
	loc::Array{Array{Float64,2},2}			#measurements
	frame::Array{Array{Int64,1},2}			#frame number corresponding to measurements
	mol_ID::Array{Array{Int64,1},2}			#population ID for each measurement
	MU::Array{Array{Float64,2},2}			#mean of measurements for each ID
	K::Array{Int64,2}	#number of molecules in segment (according to ground_truth)
	N::Array{Int64,2}	#number of observations in segment
	bounds::Array{BoundsSegment,2}	#typeof(BoundsSegment)=BoundsSegment

		#Overall (unrestricted) features
	all_ground_truth::Array{Float64}
	all_loc::Array{Float64}
	all_frame::Array{Int64}
	all_mol_ID::Array{Int64}
##Function to initialize structure
	function SegmentSR(NN::Int64,df1,df2)
		DataSim=df1
		RealLoc=df2
		obj=new()
		obj.idx=Array{Tuple{Int64,Int64},1}(undef,0)	#valid indices from restriction; use as A[CartesianIndex(obj.idx[j])]
##Segment-restricted features
		obj.ground_truth=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	ℝsim
		obj.loc=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	𝕐sim
		obj.frame=Array{Array{Int64,1},2}(undef,NN-1,NN-1)	#	𝔽sim
		obj.mol_ID=Array{Array{Int64,1},2}(undef,NN-1,NN-1)	#	𝕄sim
		obj.MU=Array{Array{Float64,2},2}(undef,NN-1,NN-1)	#	mean of members for each mol_ID
		obj.K=Array{Int64,2}(undef,NN-1,NN-1)
		obj.N=Array{Int64,2}(undef,NN-1,NN-1)
		obj.bounds=Array{BoundsSegment,2}(undef,NN-1,NN-1)
##Overall (unrestricted) features
		obj.all_ground_truth=Array{Float64}(RealLoc[[2,3,1],])	#Real Localisations
		obj.all_loc=Matrix{Float64}(DataSim[[4,5,3],])		#Observations
		obj.all_frame=Array{Int64}(DataSim[3,])	#Framenumber
		obj.all_mol_ID=Array{Int64}(DataSim[1,])	#mol_ID
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
		for j=1:(NN-1)
			for k=1:(NN-1)
				obj.bounds[j,k]=BoundsSegment([Xpart[j],Xpart[j+1]],[Ypart[k],Ypart[k+1]])
			end
		end
##3.. Restrict Area
		@showprogress for partN=1:(NN-1)
			for partN2=1:(NN-1)
				# print("VVVVVVVVVVVVVVVVVVV\n")
				# print("partN = $(partN)\n")
				# print("partN2 = $(partN2)\n")
				xmin=Xpart[partN]
				xmax=Xpart[partN+1]
				ymin=Ypart[partN2]
				ymax=Ypart[partN2+1]
				#Crop sim-data based on partition
				obj.loc[partN,partN2],val=dim2crop(𝕐sim,xmin,xmax,ymin,ymax)	#val is the boolean for indices in contraint area
				obj.mol_ID[partN,partN2]=𝕄sim[val]
				obj.frame[partN,partN2]= 𝔽sim[val]
				obj.ground_truth[partN,partN2],val2= dim2crop(ℝsim,xmin,xmax,ymin,ymax)
				push!(amountM,size(obj.loc[partN,partN2])[1])#######################substitute for object field obj.N???#########################################################################################################################
				obj.N[partN,partN2]=size(obj.loc[partN,partN2])[1]
				obj.K[partN,partN2]=size(unique(𝕄sim[val]),1)
##4..Verify if segments are valid under constraints
				if size(obj.loc[partN,partN2])[1] > 5 && size(unique(obj.mol_ID[partN,partN2]),1)>2
					push!(obj.idx,(partN,partN2))
					#Crop Framenumber, mol_ID and Real Localisations based on cropped sim-data
					#For validation purpose, determine amount of simulated molecules
					𝛍sim= Matrix{Float64}(undef,length(unique(obj.mol_ID[partN,partN2])),2)
					#For unique molecule numbers compute means
					for i=1:length(unique(obj.mol_ID[partN,partN2]))
						𝛍sim[i,:]=mean(obj.loc[partN,partN2][obj.mol_ID[partN,partN2] .== unique(obj.mol_ID[partN,partN2])[i],1:2]',dims=2)[:]
					end
					obj.MU[partN,partN2]=𝛍sim
				end
			end
		end
		return obj
	end
end
export SegmentSR
##


#
# ##
#
