using CSV
using DataFrames

function importSR(SIMurl::String,LOCurl::String)
#1.. Import Data
	#Load Data with colnames:
	#mol_ID|loc_ID|framenumber|pos_x|pos_y
	#[Int64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,Float64]
	#DataSim = CSV.read("/Users/Patrick/Documents/GitHub/BayesianEstimation/STORM/Simulationen_01-24-19/sim_Tet.csv",
	#DataSim = CSV.read("/Users/Patrick/Documents/GitHub/BayesianEstimation/STORM/ORI_SIM/locs_gamma_conf.csv",
	DataSim = CSV.read(SIMurl,
					types=[Int64,Int64,Int64,Float64,Float64],
					header=["label_ID","oligomer_ID","framenum","x_nm","y_nm"], datarow=2, DataFrame)
	#RealLoc = CSV.read("/Users/Patrick/Documents/GitHub/BayesianEstimation/STORM/Simulationen_01-24-19/sim_Tet_locs.csv",
	#RealLoc = CSV.read("/Users/Patrick/Documents/GitHub/BayesianEstimation/STORM/ORI_SIM/labels_gamma.csv",
	RealLoc = CSV.read(LOCurl,
					types=[Int64,Float64,Float64],
					header=["label_ID","x_nm","y_nm"], datarow=2, DataFrame)
	return DataSim, RealLoc
end

export importSR
