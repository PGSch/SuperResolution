using SuperResolution
using Test

@testset "segmentSR.jl" begin
    @test fieldnames(SuperResolution.SegmentSR)==(:select, :idx, :ground_truth, :loc, :frame, :mol_ID, :MU, :K, :N, :bounds, :all_ground_truth, :all_loc, :all_frame, :all_mol_ID)
    # #my_f(x,y)=2x+3y
    # @test my_f(2,1)==7
    # @test my_f(2,3)==13
    # # Write your tests here.
end
