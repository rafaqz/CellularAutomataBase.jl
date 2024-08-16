using DynamicGrids, StaticArrays, Test, FileIO, Colors, FixedPointNumbers
using DynamicGrids: SimData, NoDisplayImageOutput

@testset "CellRule that multiples a StaticArray" begin
    rule = Cell{:grid1}() do data, state, I
         2state
    end
    init = (grid1 = fill(SA[1.0, 2.0], 5, 5),)
    output = ArrayOutput(init; tspan=1:3)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0],
        SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0],
        SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], 
        SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], 
        SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0], SA[2.0, 4.0]
    ], (5, 5))
    @test output[3][:grid1] == reshape([ # Have to use reshape to construct this
        SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0],
        SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0],
        SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], 
        SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], 
        SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0], SA[4.0, 8.0]
    ], (5, 5))
end

@testset "Neighborhood Rule that sums a Neighborhood of StaticArrays" begin
    rule = Neighbors{Tuple{:grid1,:grid2},:grid1}(Moore(1)) do data, hood, (s1, s2), I
        sum(hood) .+ s2
    end
    init = (
        grid1 = fill(SA[1.0, 2.0], 5, 5),
        grid2 = fill(0.5, 5, 5),
    )
    output = ArrayOutput(init; tspan=1:2)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        SA[3.5, 6.5],  SA[5.5, 10.5], SA[5.5, 10.5], SA[5.5, 10.5], SA[3.5, 6.5],
        SA[5.5, 10.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[5.5, 10.5],
        SA[5.5, 10.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[5.5, 10.5], 
        SA[5.5, 10.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[8.5, 16.5], SA[5.5, 10.5], 
        SA[3.5, 6.5],  SA[5.5, 10.5], SA[5.5, 10.5], SA[5.5, 10.5], SA[3.5, 6.5]
    ], (5, 5))
end

@testset "SetCell randomly updates a StaticArray" begin
    rule = SetCell{:grid1}() do data, state, I
        if I == (2, 2) || I == (1, 3)
            data[:grid1][I...] = SA[99.0, 100.0]
        end
    end
    init = (grid1 = fill(SA[0.0, 0.0], 3, 3),)
    output = ArrayOutput(init; tspan=1:2)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        SA[0.0, 0.0], SA[0.0, 0.0], SA[0.0, 0.0], 
        SA[0.0, 0.0], SA[99.0, 100.0], SA[0.0, 0.0], 
        SA[99.0, 100.0], SA[0.0, 0.0], SA[0.0, 0.0]
    ], (3, 3))
end


struct TestStruct{A,B}
    a::A
    b::B
end
const TS = TestStruct

Base.:*(ts::TestStruct, x::Number) = TestStruct(ts.a * x, ts.b * x) 
Base.:*(x::Number, ts::TestStruct) = TestStruct(x * ts.a, x * ts.b)
Base.:/(ts::TestStruct, x::Number) = TestStruct(ts.a / x, ts.b / x) 
Base.:+(ts1::TestStruct, ts2::TestStruct) = TestStruct(ts1.a + ts2.a, ts1.b + ts2.b)
Base.:-(ts1::TestStruct, ts2::TestStruct) = TestStruct(ts1.a - ts2.a, ts1.b - ts2.b)
Base.:+(ts::TestStruct, x::Number) = TestStruct(ts.a + x, ts.b + x)
Base.:-(ts::TestStruct, x::Number) = TestStruct(ts.a - x, ts.b - x)

Base.isless(a::TestStruct, b::TestStruct) = isless(a.a, b.a)
Base.zero(::Type{<:TestStruct{T1,T2}}) where {T1,T2} = TestStruct(zero(T1), zero(T2))
Base.oneunit(::Type{<:TestStruct{T1,T2}}) where {T1,T2} = TestStruct(oneunit(T1), oneunit(T2))

DynamicGrids.to_rgb(scheme::ObjectScheme, obj::TestStruct) = ARGB32(obj.a)
DynamicGrids.to_rgb(scheme, obj::TestStruct) = get(scheme, obj.a)


@testset "CellRule that multiples a struct" begin
    rule = Cell{:grid1,:grid1}() do data, state, I
         2state
    end
    init = (grid1 = fill(TS(1.0, 2.0), 5, 5),)
    output = ArrayOutput(init; tspan=1:3)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0),
        TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0),
        TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), 
        TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), 
        TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0), TS(2.0, 4.0)
    ], (5, 5))
    @test output[3][:grid1] == reshape([ # Have to use reshape to construct this
        TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0),
        TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0),
        TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), 
        TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), 
        TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0)
    ], (5, 5))
end

@testset "Neighborhood Rule that sums a Neighborhood of stucts" begin
    rule = Neighbors{Tuple{:grid1,:grid2},:grid1}(Moore(1)) do data, hood, (s1, s2), I
        sum(hood) * s2
    end
    init = (
        grid1 = fill(TS(1.0, 2.0), 5, 5),
        grid2 = fill(0.5, 5, 5),
    )
    output = ArrayOutput(init; tspan=1:2)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        TS(1.5, 3.0), TS(2.5, 5.0), TS(2.5, 5.0), TS(2.5, 5.0), TS(1.5, 3.0),
        TS(2.5, 5.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(2.5, 5.0),
        TS(2.5, 5.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(2.5, 5.0), 
        TS(2.5, 5.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(4.0, 8.0), TS(2.5, 5.0), 
        TS(1.5, 3.0), TS(2.5, 5.0), TS(2.5, 5.0), TS(2.5, 5.0), TS(1.5, 3.0),
    ], (5, 5))
end

@testset "SetCell rule randomly updates a struct" begin
    rule = SetCell{:grid1,:grid1}() do data, state, I
        if I == (2, 2) || I == (1, 3)
            data[:grid1][I...] = TS(99.0, 100.0)
        end
    end
    init = (grid1 = fill(TS(0.0, 0.0), 3, 3),)
    output = ArrayOutput(init; tspan=1:2)
    sim!(output, rule)
    @test output[2][:grid1] == reshape([ # Have to use reshape to construct this
        TS(0.0, 0.0), TS(0.0, 0.0), TS(0.0, 0.0), 
        TS(0.0, 0.0), TS(99.0, 100.0), TS(0.0, 0.0), 
        TS(99.0, 100.0), TS(0.0, 0.0), TS(0.0, 0.0)
    ], (3, 3))
end

@testset "object grid can generate an image" begin
    @testset "normalise" begin
        @test DynamicGrids.to_rgb(ObjectScheme(), TestStruct(99.0, 1.0) / 99) == ARGB32(1.0)
        @test DynamicGrids.to_rgb(ObjectScheme(), TestStruct(00.0, 0.0) / 99) == ARGB32(0.0)
        @test DynamicGrids.to_rgb(ObjectScheme(), DynamicGrids.normalise(TestStruct(99.0, 1.0), nothing, 99)) == ARGB32(1.0)
    end
    rule = SetCell{:grid1,:grid1}() do data, state, I
        if I == (2, 2) || I == (1, 3)
            data[:grid1][I...] = TS(99.0, 100.0)
        end
    end
    init = (grid1=fill(TS(0.0, 0.0), 3, 3),)
    # These should have the same answer
    output1 = GifOutput(init; 
        filename="objectgrid.gif", store=true, tspan=1:2, maxval=[99.0], text=nothing
    )
    output2 = GifOutput(init; 
        filename="objectgrid_greyscale.gif", scheme=Greyscale(), store=true, tspan=1:2, 
        maxval=reshape([99.0], 1, 1), text=nothing
    )
    sim!(output1, rule)
    sim!(output2, rule)
    @test output1[2][:grid1] == 
       [TS(0.0, 0.0) TS(0.0, 0.0) TS(99.0, 100.0)
        TS(0.0, 0.0) TS(99.0, 100.0) TS(0.0, 0.0) 
        TS(0.0, 0.0) TS(0.0, 0.0) TS(0.0, 0.0)]
    @test RGB.(output1.gif[:, :, 2]) == 
          RGB.(output2.gif[:, :, 2]) == 
          load("objectgrid.gif")[:, :, 2] == 
          load("objectgrid_greyscale.gif")[:, :, 2] == 
          map(xs -> RGB{N0f8}(xs...), 
              [(0.298,0.298,0.298) (0.298,0.298,0.298) (1.0,1.0,1.0)
               (0.298,0.298,0.298) (1.0,1.0,1.0) (0.298,0.298,0.298)          
               (0.298,0.298,0.298) (0.298,0.298,0.298) (0.298,0.298,0.298)]
          )
end

@testset "static arrays grid can generate an image" begin
    rule = Cell{:grid1}() do data, state, I
         2state
    end
    init = (grid1 = fill(SA[1.0, 2.0], 5, 5),)

    # We can index into the SArray or access
    # it with a function, defined using a Pair
    output = GifOutput(init; 
        filename="sa.gif",
        tspan=1:3, 
        store=true,
        layout=[:grid1=>1 :grid1=>x->x[2]],
        scheme=[Greyscale() Greyscale()],
        minval=[0.0 0.0], maxval=[10.0 10.0],
        text=nothing,
    )

    sim!(output, rule)
    a02 = ARGB32(0.2)
    a04 = ARGB32(0.4)
    @test output.gif[:, :, 2] ==
        [a02 a02 a02 a02 a02 a04 a04 a04 a04 a04
         a02 a02 a02 a02 a02 a04 a04 a04 a04 a04
         a02 a02 a02 a02 a02 a04 a04 a04 a04 a04
         a02 a02 a02 a02 a02 a04 a04 a04 a04 a04
         a02 a02 a02 a02 a02 a04 a04 a04 a04 a04]
end

@testset "objects must be isbits" begin
    @test_throws ArgumentError ArrayOutput((grid1 = fill([1.0, 2.0], 5, 5),); tspan=1:3)
    @test_throws ArgumentError ArrayOutput((grid1=fill(1, 5, 5), grid2=fill([1.0, 2.0], 5, 5),); tspan=1:3)
    @test_throws ArgumentError ArrayOutput(fill(Ref(1), 5, 5); tspan=1:3)
    @test_throws ArgumentError TransformedOutput(identity, fill(Ref(1), 5, 5); tspan=1:3)
end
