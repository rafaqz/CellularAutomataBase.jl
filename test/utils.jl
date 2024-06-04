using DynamicGrids, Test
using DynamicGrids: inbounds, isinbounds, _inbounds, _isinbounds,
    SimData, _unwrap, ismasked

@testset "boundary boundary checks are working" begin
    @testset "inbounds with Remove() returns index and false for an boundaryed index" begin
        @test _inbounds(Remove(), (4, 5), 1, 1) == ((1,1), true)
        @test _inbounds(Remove(), (4, 5), 2, 3) == ((2,3), true)
        @test _inbounds(Remove(), (4, 5), 4, 5) == ((4,5), true)
        @test _inbounds(Remove(), (4, 5), 0, 0) == ((0,0), false)
        @test _inbounds(Remove(), (3, 2), 2, 3) == ((2,3), false)
        @test _inbounds(Remove(), (1, 4), 2, 3) == ((2,3), false)
        @test _inbounds(Remove(), (2, 3), 200, 300) == ((200,300), false)
        @test _inbounds(Remove(), (4, 5), -3, -100) == ((-3,-100), false)
    end
    @testset "inbounds with Wrap() returns new index and true for an boundaryed index" begin
        @test _inbounds(Wrap(), (10, 10),  -2, 3) == ((8,  3), true)
        @test _inbounds(Wrap(), (10, 10),   2, 0) == ((2, 10), true)
        @test _inbounds(Wrap(), (10, 10),  22, 0) == ((2, 10), true)
        @test _inbounds(Wrap(), (10, 10), -22, 0) == ((8, 10), true)
    end
    @testset "inbounds with Reflect() returns new index and true for an out of bounds index" begin
        @test _inbounds(Reflect(), (10, 10),  -2, 3) == ((4, 3), true)
        @test _inbounds(Reflect(), (10, 10),   2, 0) == ((2, 2), true)
        @test _inbounds(Reflect(), (10, 10),  12, 0) == ((8, 2), true)
        @test _inbounds(Reflect(), (10, 10), -22, 0) == ((24, 2), true)
    end
    @testset "isinbounds" begin
        @test _isinbounds((4, 5), 4, 5) == true
        @test _isinbounds((2, 3), 200, 300) == false
        @test _isinbounds((10, 10), -22, 0) == false
    end
    @testset "boundscheck objects" begin
        output = ArrayOutput(zeros(Int, 10, 10); tspan=1:10)
        sd = SimData(output.extent, Ruleset())
        @test inbounds(sd, 5, 5) == ((5, 5), true)
        @test inbounds(first(sd),  5, 5) == ((5, 5), true)
        @test inbounds(sd, 12, 5) == ((12, 5), false)
        sd_wrap = SimData(output.extent, Ruleset(; boundary=Wrap()))
        @test inbounds(sd_wrap, 5, 5) == ((5, 5), true)
        @test inbounds(sd_wrap, 12, 5) == ((2, 5), true)
        @test inbounds(first(sd_wrap), 12, 5) == ((2, 5), true)
        sd_reflect = SimData(output.extent, Ruleset(; boundary=Reflect()))
        @test inbounds(sd_reflect, 5, 5) == ((5, 5), true)
        @test inbounds(sd_reflect, 12, 5) == ((8, 5), true)
        @test inbounds(first(sd_reflect), 12, 5) == ((8, 5), true)
    end
end

@testset "isinferred" begin
    @testset "unstable conditional" begin
        rule = let threshold = 20
            Cell() do data, x, I
                x > 1 ? 2 : 0.0
            end
        end
        output = ArrayOutput(rand(Int, 10, 10); tspan=1:10)
        @test_throws ErrorException isinferred(output, rule)
    end

    @testset "return type" begin
        rule = Neighbors{:a,:a}(Moore{1}()) do data, hood, x, I
            Int(maximum(hood))
        end
        output = ArrayOutput((a=rand(Int, 10, 10),); tspan=1:10)
        @test isinferred(output, rule)
        output = ArrayOutput((a=rand(Bool, 10, 10),); tspan=1:10)
        @test_throws ErrorException isinferred(output, rule)
    end

end

@testset "ismasked" begin
    output = ArrayOutput(zeros(2, 2); mask=Bool[1 0; 0 1], tspan=1:10)
    sd = SimData(output.extent, Ruleset())
    @test ismasked(sd, 1, 2) == true
    @test ismasked(sd, 2, 2) == false
    output_nomask = ArrayOutput(zeros(2, 2); tspan=1:10)
    sd = SimData(output_nomask.extent, Ruleset())
    @test ismasked(sd, 1, 2) == false
end

@testset "unwrap" begin
    @test _unwrap(1) == 1
    @test _unwrap(Val(:a)) == :a
    @test _unwrap(Aux(:a)) == :a
    @test _unwrap(Grid(:a)) == :a
    @test _unwrap(Aux{:x}) == :x
    @test _unwrap(Grid{:x}) == :x
    @test _unwrap(Val{:x}) == :x
end
