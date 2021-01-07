using SuperResolution
using Test
using DataFrames

@testset "importSR.jl" begin
    @test typeof(SuperResolution.DataSim)==DataFrame
    @test typeof(SuperResolution.RealLocs)==DataFrame
    # Load Data with colnames:
    # mol_ID|loc_ID|framenumber|pos_x|pos_y
    # [Int64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,Float64]
    # Write your tests here.
end
