using SuperResolution
using Test

@testset "importSR.jl" begin
    @test typeof(importSR("/Users/Patrick/Documents/GitHub/SuperResolution/src/locs_gamma_conf.csv","/Users/Patrick/Documents/GitHub/SuperResolution/src/labels_gamma.csv"))==Tuple{DataFrame,DataFrame}
    # Load Data with colnames:
    # mol_ID|loc_ID|framenumber|pos_x|pos_y
    # [Int64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,Float64]
    # Write your tests here.
end
