using SuperResolution
using Test

@testset "fctSR.jl" begin
    @test true
# # function dim(A)
# # 	d1=length(A[:,1])
# # 	d2=length(A[1,:])
# # 	return d1,d2
# # end
#     @test dim(zeros(Int64,3,3))==(3,3)
#     @test dim(zeros(Int64,3,3))[1]==3
#     @test dim(zeros(Int64,3,3))[2]==3
# # function eye(d::Int64)
# # 	E=zeros(Float64,d,d)
# # 	map(i->E[i,i]=1,1:dim(E)[1])
# # 	return E
# # end
#     @test eye(3)*eye(3)==eye(3)
# # function ctrans(X::Array{Float64,2})
# # 	convert(Array{Float64,2},transpose(X))
# # end
# # function ctrans(X::Array{Array{Float64,2},1})
# # 	convert(Array{Array,1},map(i->ctrans(𝛍[1][i])[:],collect(1:length(𝛍[1]))))
# # end
#     @test ctrans([[1 2];[3 4]])==[[1 3];[2 4]]
# # function dim2clean(MU_star::Array{Array,1})
# # 	MU_star = dim2adj(MU_star)
# # 	N=maximum(size(MU_star))
# # 	𝛅=Int64(3)
# # 	for i=1:N
# # 		if isnan(MU_star[i])
# # 			for j=1:𝛅
# # 				MU_star[i,j]=-10.0
# # 			end
# # 		end
# # 	end
# # 	MU_star = dim2adj(MU_star)
# # 	MU_star = MU_star[MU_star .> dim2adj([100.0,100.0,100.0])]
# # end
# # function dim2clean(MU_star::Array{Array,1},𝛆::Array{Float64,1})
# # 	MU_star = dim2adj(MU_star)
# # 	N=maximum(size(MU_star))
# # 	𝛅=Int64(3)
# # 	for i=1:N
# # 		if isnan(MU_star[i])
# # 			for j=1:𝛅
# # 				MU_star[i,j]=-10.0
# # 			end
# # 		end
# # 	end
# # 	MU_star = dim2adj(MU_star)
# # 	MU_star = MU_star[MU_star .>= dim2adj(𝛆)]
# # end
# # function dim2crop(X::Array{Float64,2},xmin::Float64,xmax::Float64,ymin::Float64,ymax::Float64)
# # 	Xog = X
# # 	X=X[X[:,1].>=xmin,:]
# # 	X=X[X[:,1].<=xmax,:]
# # 	X=X[X[:,2].>=ymin,:]
# # 	X=X[X[:,2].<=ymax,:]
# # 	val=(Xog[:,1] .∈ [X[:,1]]) .& (Xog[:,2] .∈ [X[:,2]])
# # 	return[X,val]
# # end
# #
# # function dim2adj(X::Array{Float64,2})
# # 	J = size(X,1)
# # 	# if size(X,1)<size(X,2)
# # 	# 	J = size(X,2)
# # 	# 	X = X'
# # 	# end
# # 	MU_0 = Array{Array,1}(undef,J)
# # 	for j=1:J
# # 		MU_0[j] = X[j,:]
# # 	end
# # 	return MU_0
# # end
# # function dim2adj(X::Array{Float64,1})
# # 	MU_0 = Array{Array,1}(undef,1)
# # 	MU_0[1] = X[:,1]
# # 	return MU_0
# # end
# # function dim2adj(X::Array{Array,1},δ=Int64(2))
# # 	if δ != Int64(2)
# # 		Y_raw = Array{Float64,2}(undef,maximum(size(X)),size(X[1],1))
# # 		for i=1:maximum(size(X))
# # 			for j=1:size(X[1],1)
# # 				Y_raw[i,j] = X[i][j]
# # 			end
# # 		end
# # 		return Y_raw
# # 	end
# # 	Y_raw = Array{Float64,2}(undef,maximum(size(X)),2)
# # 	for i=1:maximum(size(X))
# # 		Y_raw[i,1] = X[i][1]
# # 		Y_raw[i,2] = X[i][2]
# # 	end
# # 	return Y_raw
# # end
# #
# # ####Generate MixtureModel off given data
# # function dim2mix(X::Array{Array,1},S=eye(size(X,2)))
# # 	MM = MixtureModel(map(x->MvNormal(x,S),X))
# # 	return MM
# # end
# # function dim2mix(X::Array{Float64,2},S=eye(size(X,2)))
# # 	X = dim2adj(X)
# # 	MM = MixtureModel(map(x->MvNormal(x,S),X))
# # 	return MM
# # end
# # function dim2mix(X::Array{Float64,2},w::Array{Float64,1},S=eye(size(X,2)))
# # 	#print("[dim2mix] X=",X,"\n")
# # 	X = dim2adj(X)
# # 	#print("[dim2mix..dim2adj(X)] X=",X,"\n")
# # 	#print("[dim2mix] S=",S,"\n")
# # 	MM = MixtureModel(map(x->MvNormal(x,S),X))
# # 	return MM
# # end
# # function dim2mix(X::Array{Array,1},Σ::Array{Array{Float64,2}})
# # 	MM = MixtureModel(map(i->MvNormal(X[i],Σ[i]),1:length(X)))
# # 	return MM
# # end
# # function dim2mix(X::Array{Array,1},Σ::Array{Float64,2},p::Array{Float64,1})
# # 	MM = MixtureModel(map(x->MvNormal(x,Σ),X),p)
# # 	return MM
# # end
# # function dim2mix(X::Array{Float64,2},Σ::Array{Float64,2})
# # 	X = dim2adj(X)
# # 	MM = MixtureModel(map(x->MvNormal(x,Σ),X))
# # 	return MM
# # end
# # ####Evaluate bounds of data and generate approriate range
# # ####in bounds for package KernelDensity
# # function dim2bounds(X::Array{Float64,2})
# # 	M = 0*Matrix{Float64}(undef,2,2)
# # 	M[1,1] = minimum(X[:,1])
# # 	M[2,1] = maximum(X[:,1])
# # 	M[1,2] = minimum(X[:,2])
# # 	M[2,2] = maximum(X[:,2])
# # 	return M
# # end
# # function dim2range(strt::Float64,stp::Float64,len::Int64)
# # 	Omega = range(strt,stop=stp,length=Int64(len))
# # 	return Omega
# # end
# # ####Generate appropriate range in bounds for data X_raw in dimension d
# # ####The range has to be generated for each dimension
# # function dim2range(X_raw::Array{Float64,2},d::Int64)
# # 	if size(X_raw,1)<size(X_raw,2)
# # 		X_raw=ctrans(X_raw)
# # 	end
# # 	len = Int64(floor(sqrt(maximum(size(X_raw)))))
# # 	M = dim2bounds(X_raw)
# # 	Omega = range(M[1,d],stop=M[2,d],length=Int64(len))
# # 	return Omega
# # end
# # function dim2range(X_raw::Array{Float64,2})
# # 	if size(X_raw,1)<size(X_raw,2)
# # 		X_raw=ctrans(X_raw)
# # 	end
# # 	len = Int64(floor(sqrt(maximum(size(X_raw)))))
# # 	M = dim2bounds(X_raw)
# # 	Omega1 = range(M[1,1],stop=M[2,1],length=Int64(len))
# # 	Omega2 = range(M[1,2],stop=M[2,2],length=Int64(len))
# # 	delta = [Omega1[end]-Omega1[1];Omega2[end]-Omega2[1]]
# # 	return [Omega1,Omega2,delta]
# # end
# # #MLE geom dis given FrNm
# # function geomMLE(FrNm::Array{Int64,1})
# # 	β=collect(1e-8:1e-2:1.0)
# # 	βpost=[]
# # 	#Likelihood for β parameter of Geometric Distribution, based on Framenmbers of Observations
# # 	FrNm = FrRef(FrNm)
# # 	for j=1:length(β)
# # 		push!(βpost,sum([FrNm[i]*log(1-β[j])+log(β[j])] for i in 1:length(FrNm))[1])
# # 	end
# # 	#MLE for β Parameter
# # 	βhat=β[collect(1:length(β))[βpost.==maximum(βpost)]][1]
# # 	#Geometric Distribution with MLE for Parameter
# # 	G=Geometric(βhat)
# # 	return G
# # end
# # function normalMLE(FrNm::Array{Int64,1})
# # 	μMLE=1/length(FrNm) * sum(FrNm)
# # 	σ=1.0
# # 	nMLE=Normal(μMLE,σ)
# # 	# σMLE=sqrt(1/length(FrNm) * sum( (FrNm .- μMLE).^2 ) )
# # 	# if sum(FrNm .== μMLE) == length(FrNm) σMLE = 1.0 end #avoid 0 value for σMLE (occurs for all FrNm==μMLE)
# # 	#nMLE=Normal(μMLE,σMLE)
# # 	return nMLE
# # end
# # ##ref FrNm of obs to first frame
# # function FrRef(FrNm::Array{Int64,1})
# # 	FrNm = sort(FrNm) .- minimum(FrNm)
# # end
# # #
# # function FrRef(FrNm::Int64)
# # 	return FrNm
# # end
# # #compute data based on cluster corresponding assignment in AcollElem
# # function Feval(AcollElem::Array{Int64,1})#,𝔽sim::Array{Int64,1})
# # 	uniqueLengthElem=length(unique(AcollElem))
# # 	uniqueAcollElem=unique(AcollElem)
# # 	#get framenumbers for each cluster
# # 	𝔽evalElem=Array{Array,1}(undef,uniqueLengthElem)
# # 	𝔽evalElem=map(i->𝔽sim[AcollElem .== unique(AcollElem)[i]],1:uniqueLengthElem)
# # 	#eval indicies of observations (to evaluate locations in 𝕐sim)
# # 	𝔽evalElemIdx=Array{Array,1}(undef,uniqueLengthElem)
# # 	𝔽evalElemIdx=map(i->findall(x->x∈𝔽evalElem[i],𝔽sim),1:uniqueLengthElem)
# # 	#get observations corresponding to clusters
# # 	𝕐evalElem=Array{Array,1}(undef,uniqueLengthElem)
# # 	𝕐evalElem=map(i->𝕐sim[𝔽evalElemIdx[i],:],1:uniqueLengthElem)
# # 	#FrmNm f.e. Cl ID||Indices of obs f.e. Cl ID||Loc of obs f.e. Cl ID||number of unique Cl||ID of unique Cl
# # 	return 𝔽evalElem,𝔽evalElemIdx,𝕐evalElem,uniqueLengthElem,uniqueAcollElem
# # end
# # #Kullback Leibler
# # function mu_star(mu1::Array{Float64,1},mu2::Array{Float64,1},W::Array{Float64,1})
# # 	mu_star = (W[1]*mu1+W[2]*mu2)/(W[1]+W[2])
# # 	return mu_star
# # end
# # function s_star(mu1::Array{Float64,1},mu2::Array{Float64,1},S1::Array{Float64,2},S2::Array{Float64,2},W::Array{Float64,1})
# # 	Σ_star = (W[1]*S1 + W[2]*S2)/(W[1]+W[2]) + (W[1]*W[2]*(mu1-mu2)*(mu1-mu2)')/((W[1]+W[2])^2)
# # 	return Σ_star
# # end
# # function dim2KLD(mu1::Array{Float64,1},mu2::Array{Float64,1},S1::Array{Float64,2},S2::Array{Float64,2},W::Array{Float64,1})
# # 	S=[S1,S2]
# # 	Σ_star = s_star(mu1,mu2,S[1],S[2],[W[1],W[2]])
# # 	𝛅=Int64(3)
# # 	KLD = (sum([W[i]/2.0 * (tr(Σ_star^-1 * S[i]) + log(det(Σ_star * S[i]^-1))-𝛅)] for i in 1:2)[1] + (W[1]*W[2]/(2.0*(W[1]+W[2])^2) * (mu1-mu2)'*Σ_star^-1*(mu1-mu2)))
# # 	if isnan(KLD) KLD=0.0 end
# # 	return KLD
# # end
# # #Extrema of KLD matrix
# # function dim2KLDex(KLDfull::Array{Float64,2})
# # 	𝛕=1e-5
# # 	EX = 0*Array{Int64,2}(undef,3,2)
# # 	K = maximum(size(KLDfull))
# # 	#column1 minimum KLD indices [vectorindex, matrixindices[j,k]]
# # 	EX[2,1] = findall(KLDfull .== minimum(KLDfull[findall((KLDfull .> 𝛕))]))[1][1] #column index
# # 	EX[3,1] = findall(KLDfull .== minimum(KLDfull[findall((KLDfull .> 𝛕))]))[1][2] #row index
# # 	EX[1,1] = Int64(floor((EX[3,1]-1)*K+EX[2,1])) #single indexing
# # 	EXmin = KLDfull[EX[1,1]]
# # 	#column2 maximum KLD indices [vectorindex, matrixindices[j,k]]
# # 	EX[2,2] = findall(KLDfull .== maximum(KLDfull[findall((KLDfull .> 𝛕))]))[1][1]
# # 	EX[3,2] = findall(KLDfull .== maximum(KLDfull[findall((KLDfull .> 𝛕))]))[1][2]#Int64(EX[1,2]-K*(EX[3,2]-1))
# # 	EX[1,2] = Int64(floor((EX[3,2]-1)*K+EX[2,2]))
# # 	EXmean=mean(collect(KLDfull[:]))
# # 	#print("Idx min: ",EX[1,1],"= ","[",EX[2,1],",",EX[3,1],"]","min value: ", EXmin,"(mean value: ",EXmean,")\n")
# # 	return EX
# # end
# # ##Weighted Kullback Leibler Divergence
# # function dim2WKLDfull(MU::Array{Array,1},S::Array{Array{Float64,2},1},W::Array{Float64,1},A=NaN,Atmp=NaN)
# # 	#entries are merged during each collapse, therefore K represents the current number of clusters
# # 	K=maximum(size(MU))
# # 	#if no A or Atmp is submitted, it's set as NaN (though A and Atmp are by default submitted during dim2WKLDfull computation)
# # 	if isnan(A[1]) A=collect(1:K) end
# # 	if isnan(Atmp[1]) Atmp=collect(1:K) end
# # 	KLD = 1e3*ones(K,K)
# # 	global MU_star = map(i->MU[i],1:length(MU)) #(clustered-)data based on 𝕐sim
# # 	global A_star = map(i->A[i],1:length(A)) #cluster ID
# # 	#S_star = S
# # 	W_star = map(i->W[i],1:length(W)) #weights
# # 	A_coll = map(i->Atmp[i],1:length(Atmp)) #A_coll stores the cluster IDs which are updated for every collapse
# # 	𝛅 = Int64(3) #dimension
# # 	S_star= repeat([1.0*eye(𝛅)],length(S)) #COV Matrices
# # 	FrNm = Feval(A_coll)[1] #evaluate frames f.e. clust ID in A_coll
# # 	#compute normal distribution with μMLE and σMLE for each current cluster
# # 	#NormDistClustID = map(i->normalMLE(FrNm[i]),1:length(FrNm))
# # 	GeomDistClustID = map(i->geomMLE(FrRef(FrNm[i])),1:length(FrNm))
# # 	#KL distances
# # 	K = maximum(size(MU_star))
# # 	#print("K=",K,"\n")
# # 	#initialize "distance" matrix
# # 	KLD = 1e-15*ones(K,K)
# # 	#compute KL-distances for each and every pair of points
# # 	for k=1:K
# # 		for j=1:K
# # 			if j==k
# # 				KLD[j,k]=1e5
# # 			else
# # 				#KLD[j,k] = (2-sum(pdf.(NormDistClustID[k],FrNm[j]))) * dim2KLD(MU_star[k],MU_star[j],S_star[k],S_star[j],[1.0,1.0])
# # 				KLD[j,k] = (1/exp(sum(pdf.(GeomDistClustID[k],FrNm[j])))) * dim2KLD(MU_star[k],MU_star[j],S_star[k],S_star[j],[1.0,1.0])
# # 			end
# # 		end
# # 	end
# # 	#dim2KLDex computes extrema and corresponding indices in the distance matrix
# # 	KLD[map(i->isnan(KLD[i]),1:length(KLD[:]))].= maximum(KLD[map(i->!isnan(KLD[i]),1:length(KLD[:]))])
# # 	EX = dim2KLDex(KLD)
# # 	#min1,min2 is the pair of values with minimal KLD and get merged
# # 	min1=EX[2,1]
# # 	min2=EX[3,1]
# # 	global MU_star[EX[2,1]] = mu_star(MU_star[EX[2,1]],MU_star[EX[3,1]],[W_star[EX[2,1]],W_star[EX[3,1]]])
# # 	S_star[EX[2,1]] = s_star(MU_star[EX[2,1]],MU_star[EX[3,1]],S_star[EX[2,1]],S_star[EX[3,1]],[W_star[EX[2,1]],W_star[EX[3,1]]])
# # 	W_star[EX[2,1]] += W_star[EX[3,1]]
# # 	S_star[EX[3,1]] = NaN.*S_star[EX[3,1]]
# # 	global MU_star[EX[3,1]] =NaN.*MU_star[EX[3,1]]
# # 	W_star[EX[3,1]] = NaN.*W_star[EX[3,1]]
# # 	A_coll[A_coll .== A_star[EX[3,1]]].= A_star[EX[2,1]]
# # 	#global A_star[EX[3,1]] = 0*A_star[EX[3,1]]
# # 	global A_star[EX[3,1]] = 0*A_star[EX[3,1]]
# # 	global A_star=A_star[map(i->A_star[i]!=0,1:length(A_star))]
# # 	S_star=S_star[map(i->!isnan(S_star[i][1]),1:length(S_star))]
# # 	global MU_star=MU_star[map(i->!isnan(MU_star[i][1]),1:length(MU_star))]
# # 	W_star=W_star[map(i->!isnan(W_star[i]),1:length(W_star))]
# # 	S_star=repeat([1.0*eye(𝛅)],length(S_star))
# # 	KLD[EX[1,1]]=0.0
# # 	#print("END OF while loop dim2WKLDfull\n")
# # 	return MU_star,W_star,S_star,A_star,A_coll,KLD
# # end
end
