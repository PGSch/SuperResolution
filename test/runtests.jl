using SuperResolution
using SafeTestsets
@safetestset "GAP Tests"                begin include("clusGAP_tests.jl") end   #gap statistic
@safetestset "SR tests"                 begin include("clusSR_tests.jl") end  #interaction of all functions
@safetestset "Evaluation Tests"         begin include("evalSR_tests.jl") end    #Diagnostics
@safetestset "Misc Func Tests"          begin include("fctSR_tests.jl") end #mostly functions to change var-type/dimensions
@safetestset "Import Tests"             begin include("importSR_tests.jl") end  #data imported correctly
@safetestset "Model Tests"              begin include("modelSR_tests.jl") end    #divergence and parameter iteration
@safetestset "Plot Tests"               begin include("plotSR_tests.jl") end    #PyPlot
@safetestset "Partitioning Tests"       begin include("segmentSR_tests.jl") end    #partitioning
