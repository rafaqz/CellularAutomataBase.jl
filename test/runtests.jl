using DynamicGrids, Aqua, SafeTestsets

Aqua.test_all(DynamicGrids; ambiguities=false)

@time @safetestset "generated" begin include("generated.jl") end
@time @safetestset "rules" begin include("rules.jl") end
@time @safetestset "simulationdata" begin include("simulationdata.jl") end
@time @safetestset "utils" begin include("utils.jl") end
@time @safetestset "chain" begin include("wrappers/chain.jl") end
@time @safetestset "condition" begin include("wrappers/condition.jl") end
@time @safetestset "combine" begin include("wrappers/combine.jl") end
@time @safetestset "outputs" begin include("outputs.jl") end
@time @safetestset "transformed" begin include("transformed.jl") end
@time @safetestset "integration" begin include("integration.jl") end
@time @safetestset "objectgrids" begin include("objectgrids.jl") end
@time @safetestset "parametersources" begin include("parametersources.jl") end
@time @safetestset "show" begin include("show.jl") end
@time @safetestset "textconfig" begin include("textconfig.jl") end
# ImageMagick breaks in windows travis for some reason
if !Sys.iswindows() 
    @time @safetestset "image" begin include("image.jl") end
    @time @safetestset "makie" begin include("makie.jl") end
end
