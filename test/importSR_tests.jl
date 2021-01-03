using SuperResolution
using Test

@testset "importSR.jl" begin
    @test my_e(1,1)==2
    # Load Data with colnames:
    # mol_ID|loc_ID|framenumber|pos_x|pos_y
    # [Int64,Int64,Int64,Float64,Float64,Float64,Float64,Float64,Float64,Float64]
    # Write your tests here.
end
