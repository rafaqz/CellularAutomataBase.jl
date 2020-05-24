using DynamicGrids, Setfield, FieldMetadata, Test
import DynamicGrids: applyrule, applyrule!, maprule!, 
       source, dest, currenttime, getdata, combinedata, ruleloop,
       SimData, WritableGridData, _Read_, _Write_,
       Rule, readkeys, writekeys, @Image, @Graphic, @Output

init  = [0 1 1 0;
         0 1 1 0;
         0 1 1 0;
         0 1 1 0;
         0 1 1 0]

@testset "Ruleset or sim must include an init grid" begin
    output = ArrayOutput(init, 7)
    ruleset = Ruleset()
    @test_throws ArgumentError sim!(output, ruleset)
end


struct AddOneRule{R,W} <: Rule{R,W} end
DynamicGrids.applyrule(::AddOneRule, data, state, args...) = state + 1

@testset "Rulset mask ignores false cells" begin
    init = [0.0 4.0 0.0
            0.0 5.0 8.0
            3.0 6.0 0.0]

    mask = Bool[0 1 0
                0 1 1
                1 1 0]

    rules = Ruleset(AddOneRule{:_default_,:_default_}(); init=init, mask=mask)

    output = ArrayOutput(init, 3)
    sim!(output, rules; tspan=(1, 3))

    @test output[1] == [0.0 4.0 0.0
                        0.0 5.0 8.0
                        3.0 6.0 0.0]
    @test output[2] == [0.0 5.0 0.0
                        0.0 6.0 9.0
                        4.0 7.0 0.0]
    @test output[3] == [0.0 6.0 0.0
                        0.0 7.0 10.0
                        5.0 8.0 0.0]
end


# Single grid rules

struct TestRule{R,W} <: Rule{R,W} end
applyrule(::TestRule, data, state, index) = 0

@testset "A rule that returns zero gives zero outputs" begin
    final = [0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0
             0 0 0 0]

    rule = TestRule{:a,:a}()
    ruleset = Ruleset(rule; init=(a=init,))

    @test DynamicGrids.init(ruleset)[:a] === init
    @test DynamicGrids.mask(ruleset) === nothing
    @test DynamicGrids.overflow(ruleset) === RemoveOverflow()
    @test DynamicGrids.opt(ruleset) === SparseOpt()
    @test DynamicGrids.cellsize(ruleset) === 1
    @test DynamicGrids.timestep(ruleset) === 1
    @test DynamicGrids.ruleset(ruleset) === ruleset

    simdata = SimData(init, ruleset, 1)

    # Test maprules components
    rkeys, rdata = getdata(_Read_(), rule, simdata)
    wkeys, wdata = getdata(_Write_(), rule, simdata)
    @test rkeys == Val{:_default_}()
    @test wkeys == Val{:_default_}()
    newsimdata = @set simdata.data = combinedata(rkeys, rdata, wkeys, wdata)
    @test newsimdata.data[1] isa WritableGridData
    # Test type stability
    @inferred ruleloop(NoOpt(), rule, newsimdata, rkeys, rdata, wkeys, wdata)
    
    resultdata = maprule!(simdata, rule)
    @test source(resultdata[:_default_]) == final
end

struct TestPartial{R,W} <: PartialRule{R,W} end
applyrule!(::TestPartial, data, state, index) = 0

@testset "A partial rule that returns zero does nothing" begin
    rule = TestPartial()
    ruleset = Ruleset(rule; init=init)
    # Test type stability
    simdata = SimData(init, ruleset, 1)
    rkeys, rdata = getdata(_Read_(), rule, simdata)
    wkeys, wdata = getdata(_Write_(), rule, simdata)
    newsimdata = @set simdata.data = combinedata(wkeys, wdata, rkeys, rdata)

    @inferred ruleloop(NoOpt(), rule, newsimdata, rkeys, rdata, wkeys, wdata)

    resultdata = maprule!(simdata, rule)
    @test source(resultdata[:_default_]) == init
end

struct TestPartialWrite{R,W} <: PartialRule{R,W} end
applyrule!(::TestPartialWrite, data, state, index) = data[:_default_][index[1], 2] = 0

@testset "A partial rule that writes to dest affects output" begin
    final = [0 0 1 0;
             0 0 1 0;
             0 0 1 0;
             0 0 1 0;
             0 0 1 0]

    rule = TestPartialWrite()
    ruleset = Ruleset(rule; init=init)
    simdata = SimData(init, ruleset, 1)
    resultdata = maprule!(simdata, rule)
    @test source(first(resultdata)) == final
end

struct TestCellTriple{R,W} <: CellRule{R,W} end
applyrule(::TestCellTriple, data, state, index) = 3state

struct TestCellSquare{R,W} <: CellRule{R,W} end
applyrule(::TestCellSquare, data, (state,), index) = state^2

@testset "Chained cell rules work" begin
    init  = [0 1 2 3;
             4 5 6 7]

    final = [0 9 36 81;
             144 225 324 441]
    rule = Chain(TestCellTriple(), 
                 TestCellSquare())
    ruleset = Ruleset(rule; init=init)
    simdata = SimData(init, ruleset, 1)
    resultdata = maprule!(simdata, ruleset.rules[1]);
    @test source(first(resultdata)) == final
end


struct PrecalcRule{R,W,P} <: Rule{R,W}
    precalc::P
end
DynamicGrids.precalcrules(rule::PrecalcRule, data) = 
    PrecalcRule(currenttime(data))
applyrule(rule::PrecalcRule, data, state, index) = rule.precalc[]

@testset "Rule precalculations work" begin
    init  = [1 1;
             1 1]

    out2  = [2 2;
             2 2]

    out3  = [3 3;
             3 3]

    rule = PrecalcRule(1)
    ruleset = Ruleset(rule; init=init)
    output = ArrayOutput(init, 3)
    sim!(output, ruleset; tspan=(1, 3))
    @test output[2] == out2
    @test output[3] == out3
end


# Multi grid rules

struct DoubleY{R,W} <: CellRule{R,W} end
applyrule(rule::DoubleY, data, (x, y), index) = y * 2

struct HalfX{R,W} <: CellRule{R,W} end
applyrule(rule::HalfX, data, x, index) = x, x * 0.5

struct Predation{R,W} <: CellRule{R,W} end
Predation(; prey=:prey, predator=:predator) = 
    Predation{Tuple{predator,prey},Tuple{prey,predator}}()
applyrule(::Predation, data, (predators, prey), index) = begin
    caught = 2predators
    # Output order is the reverse of input to test that can work
    prey - caught, predators + caught * 0.5
end

predation = Predation(; prey=:prey, predator=:predator)

@testset "Multi-grid keys are inferred" begin
    @test writekeys(predation) == (:prey, :predator)
    @test readkeys(predation) == (:predator, :prey)
    @test keys(predation) == (:prey, :predator)
    @inferred writekeys(predation)
    @inferred readkeys(predation)
    @inferred keys(predation)
end

@testset "Multi-grid keys are inferred" begin
    @test writekeys(predation) == (:prey, :predator)
    @test readkeys(predation) == (:predator, :prey)
    @test keys(predation) == (:prey, :predator)
    @inferred writekeys(predation)
    @inferred readkeys(predation)
    @inferred keys(predation)
end

@testset "Multi-grid rules work" begin
    init = (prey=[10. 10.], predator=[1. 0.])
    ruleset = Ruleset(DoubleY{Tuple{:predator,:prey},:prey}(), predation)
    output = ArrayOutput(init, 3)
    sim!(output, ruleset; init=init, tspan=(1, 3))
    @test output[2] == (prey=[18. 20.], predator=[2. 0.])
    @test output[3] == (prey=[32. 40.], predator=[4. 0.])
end

@testset "Multi-grid rules work" begin
    init = (prey=[10. 10.], predator=[0. 0.])
    ruleset = Ruleset(HalfX{:prey,Tuple{:prey,:predator}}())
    output = ArrayOutput(init, 3)
    sim!(output, ruleset; init=init, tspan=(1, 3))
    @test output[2] == (prey=[10. 10.], predator=[5. 5.])
    @test output[3] == (prey=[10. 10.], predator=[5. 5.])
end

@testset "life with generic constructors" begin
    @test Life(RadialNeighborhood{1}(), (1, 1), (5, 5)) ==
          Life(; neighborhood=RadialNeighborhood{1}(), b=(1, 1), s=(5, 5))
    @test Life{:a,:b}(RadialNeighborhood{1}(), (1, 1), (5, 5)) ==
          Life(; read=:a, write=:b, neighborhood=RadialNeighborhood{1}(), b=(1, 1), s=(5, 5));
    @test Life(read=:a, write=:b) == Life{:a,:b}()
    @test Life() == Life(; read=:_default_)
end

@testset "generic ConstructionBase compatability" begin
    life = Life{:x,:y}(; neighborhood=RadialNeighborhood{2}(), b=(1, 1), s=(2, 2))
    @set! life.b = (5, 6)

    @test life.b == (5, 6)
    @test life.s == (2, 2)
    @test readkeys(life) == :x
    @test writekeys(life) == :y
    @test DynamicGrids.neighborhood(life) == RadialNeighborhood{2}()
end