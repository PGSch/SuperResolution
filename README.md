# SuperResolution

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pgsch.github.io/SuperResolution/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pgsch.github.io/SuperResolution/)
[![Build Status](https://travis-ci.com/PGSch/SuperResolution.svg?branch=master)](https://travis-ci.com/PGSch/SuperResolution.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/lsftvomwnugjr87a?svg=true)](https://ci.appveyor.com/project/PGSch/superresolution)
[![codecov](https://codecov.io/gh/PGSch/SuperResolution/branch/master/graph/badge.svg?token=zdFsH1KZKH)](https://codecov.io/gh/PGSch/SuperResolution)

This package supports the identification of blinking artifacts in super resolution microscopy.

# 1. Overview
## 1.1 Import Data
```
filename1="locs_gamma_conf.csv"
filename2="labels_gamma.csv"
filepath1 = joinpath(@__DIR__, filename1)
filepath2 = joinpath(@__DIR__, filename2)
global DataSim, RealLocs=importSR(filepath1,filepath2)
```
Where
```
function importSR(SIMurl::String,LOCurl::String)
#1.. Import Data
	#Load Data with colnames:
	#mol_ID|loc_ID|framenumber|pos_x|pos_y
	#[Int64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,Float64]
	DataSim = CSV.read(SIMurl,
					types=[Int64,Int64,Int64,Float64,Float64],
					header=["label_ID","oligomer_ID","framenum","x_nm","y_nm"], datarow=2, DataFrame)
	RealLoc = CSV.read(LOCurl,
					types=[Int64,Float64,Float64],
					header=["label_ID","x_nm","y_nm"], datarow=2, DataFrame)
	return DataSim, RealLoc
end
```


## 1.2 Segment Data into 100x100 segments
Divide data into 100x100 (SEG_SPLIT) segments, and select 5 of those randomly (SEG_RESTR) for test purpose.
```
SEG_SPLIT=100;
tmp1=SuperResolution.SegmentSR(100,DataSim,RealLocs)

```
Where
```
mutable struct SegmentSR
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
	(...)
	(...)
end
```
## Example:
```
tmp1.ground_truth[3]
```
returns the ground truth localisation of populations in segment No.3.
```
tmp1.loc[5]
```
returns all measurements [x,y,FrNm] for segment No.5 .
```
tmp1.MU[3]
```
returns the means for each population in segment No.3 .
For test purpose it might be helpful to select a small subset of all segments:
```
NUM_SEG=5;
SEG_RESTR=rand(tmp1.idx,NUM_SEG);
```
Here we randomly select 5 out of all generated segments through SegmentSR(), and save the array of indices as
```
SEG_RESTR::Array{Tuple{Int64,Int64},1}
```


## 1.3 Run Clustering
After importing, segmenting and selecting a small subset for testing, we can finally run the clustering by calling the EvalSR() function:
```
SuperResolution.EvalSR(tmp1,SEG_RESTR)
```
where
```
function EvalSR(SEG::SuperResolution.SegmentSR,SEG_RESTR::Array{Tuple{Int64,Int64},1})
```
and
```
mutable struct EvalSR
    _selectIDX::Array{Tuple{Int64,Int64},1}
    wModel::Array{SuperResolution.ModelSR,1}
    RefModel::Array{SuperResolution.ModelSR,1}
    distSR::Array{SuperResolution.DistanceSR,1}
    gapSR::Array{SuperResolution.clusGAP,1}
    diff_WK_max_idx::Array{Int64,1}
    ground_truth_K::Array{Int64,1}
    k_corr_rate::Float64

		(...)
end
```
- `wModel` for segment create all WKLD models. Start at state where every measurement is a cluster, then merge and generate model for each step until all measurements are assigned to one single cluster.

- `RefModel` same procedure as in wModel, but here use KLD clustering without the use of a weight-parameter as a reference to the weighted model.

- `distSR` the type DistanceSR, is an array with entries of type SegDistanceSR which contains information about the dissimilarity (dissimilarity matrix) where
```
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

		(...)
end
```
- `gapSR` contains the GAP-statistics.
- `diff_WK_max_idx` is a copy of the corresponding entry in the gapSR type, giving the index of solution according to GAP-statistics.
- `ground_truth` is again just a copy from the import data
- `k_corr_rate` represents the percentage difference of estimated molecules from the ground truth:
```
k_corr_rate = abs(estimated-ground_truth_number)/ground_truth_number
```
# 2. Elaborations
## 2.1 GAP-Statistics

## 2.2 WKLD

## 2.3 Recursive Parameter Estimation
