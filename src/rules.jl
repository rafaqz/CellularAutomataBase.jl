const Neighborhood = Stencils.Stencil

"""
    Rule

A `Rule` object contains the information required to apply an
`applyrule` method to every cell of every timestep of a simulation.


The [`applyrule`](@ref) method follows the form:

```julia
@inline applyrule(data::AbstractSimData, rule::MyRule, state, I::Tuple{Int,Int}) = ...
```

Where `I` is the cell index, and `state` is a single value, or a `NamedTuple`
if multiple grids are requested. the [`AbstractSimData`](@ref) object can be used to access 
current timestep and other simulation data and metadata.

Rules can be updated from the original rule before each timestep, in [`modifyrule`](@ref).
Here a paremeter depends on the sum of a grid:

```jldoctest 
using DynamicGrids, Setfield
struct MySummedRule{R,W,T} <: CellRule{R,W}
    gridsum::T
end
function modifyrule(rule::MySummedRule{R,W}, data::AbstractSimData) where {R,W}
    Setfield.@set rule.gridsum = sum(data[R])
end

# output
modifyrule (generic function with 1 method)
```

Rules can also be run in sequence, as a `Tuple` or in a [`Ruleset`](@ref)s.

DynamicGrids guarantees that:

- `modifyrule` is run once for every rule for every timestep.
    The result is passed to `applyrule`, but not retained after that.
- `applyrule` is run once for every rule, for every cell, for every timestep, unless an
    optimisation like [`SparseOpt`](@ref) is used to skip empty cells.
- the output of running a rule for any cell does not affect the input of the
    same rule running anywhere else in the grid.
- rules later in the sequence are passed grid state updated by the earlier rules.
- masked areas, and wrapped or removed `boundary` regions are updated between rules when
    they have changed.

## Multiple grids

The keys of the init `NamedTuple` will be match the grid keys used in `R` and `W` for each
rule, which is a type like `Tuple{:key1,:key1}`. Note that the names are user-specified,
and should never be fixed by a `Rule`.

Read grid names are retrieved from the type here as `R1` and `R2`, while write grids are
`W1` and `W2`.

```jldoctest
using DynamicGrids
struct MyMultiSetRule{R,W} <: SetCellRule{R,W} end
function applyrule(
    data::AbstractSimData, rule::MyMultiSetRule{Tuple{R1,R2},Tuple{W1,W2}}, (r1, r2), I
) where {R1,R2,W1,W2}
    add!(data[W1], 1, I) 
    add!(data[W2], 1, I) 
end

# output
applyrule (generic function with 1 method)
```

The return value of an `applyrule` is written to the current cell in the specified `W`
write grid/s. `Rule`s writing to multiple grids simply return a `Tuple` in the order
specified by the `W` type params.

## Rule Performance

Rules may run many millions of times during a simulation. They need to be fast.
Some basic guidlines for writing rules are:

- Never allocate memory in a `Rule` if you can help it.
- Type stability is essential. [`isinferred`](@ref) is useful to check
    if your rule is type-stable.
- Using the `@inline` macro on `applyrule` can help force inlining your
    code into the simulation.
- Reading and writing from multiple grids is expensive due to additional load
    on fast cahce memory. Try to limit the number of grids you use.
- Use a graphical profiler, like ProfileView.jl, to check your rules overall
    performance when run with `sim!`.
"""
abstract type Rule{R,W} end

# Rules are also traits - because we use wrapper rules
ruletype(::Rule) = Rule

# Allow specifying CartesianIndex for `I`
@inline applyrule(data, rule, state, I) = applyrule(data, rule, state, CartesianIndex(I))
# # Allow dropping the `data` argument
@inline applyrule(data, rule, state, I::CartesianIndex) = applyrule(rule, state, Tuple(I))
# # Allow dropping the `data` argument and using CartesianIndex
@inline applyrule(rule, state, I) = applyrule(rule, state, CartesianIndex(I))
# # Allow dropping the `I` argument
@inline applyrule(rule, state, I::CartesianIndex) = applyrule(rule, state)
# TODO add a specific method that includes all the above possible signiatures
applyrule(rule::R, state::T) where {R,T} = throw(MethodError(applyrule, Tuple{R,T}))

#= Default constructors for all Rules.
Sets both the read and write grids to `:_default`.

This strategy relies on a one-to-one relationship between fields
and type parameters, besides the initial `R` and `W` params.  =#

# No {R,W} with args or kw
function (::Type{T})(args...; kw...) where T<:Rule
    T{DEFAULT_KEY,DEFAULT_KEY}(args...; kw...)
end
# {R,W} with args
function (::Type{T})(args...) where T<:Rule{R,W} where {R,W}
    # _checkfields(T, args)
    T{map(typeof, args)...}(args...)
end
# Only R specified
function (::Type{T})(args...; kw...) where T<:Rule{R} where R
    T{R}(args...; kw...)
end

@generated function Base.keys(rule::Rule{R,W}) where {R,W}
    Expr(:tuple, QuoteNode.(union(_asiterable(W), _asiterable(R)))...)
end

@inline _writekeys(rule::Rule) = _writekeys(typeof(rule))
@inline _writekeys(::Type{<:Rule{R,W}}) where {R,W} = W
@generated function _writekeys(::Type{<:Rule{R,W}}) where {R,W<:Tuple}
    Expr(:tuple, QuoteNode.(W.parameters)...)
end

@inline _readkeys(rule::Rule) = _readkeys(typeof(rule))
@inline _readkeys(::Type{<:Rule{R,W}}) where {R,W} = R
@generated function _readkeys(::Type{<:Rule{R,W}}) where {R<:Tuple,W}
    Expr(:tuple, QuoteNode.(R.parameters)...)
end

# Define the constructor for generic rule reconstruction in Flatten.jl and Setfield.jl
function ConstructionBase.constructorof(::Type{T}) where T<:Rule{R,W} where {R,W}
    T.name.wrapper{R,W}
end

# Find the largest radius present in the passed in rules.
function radius(rules::Tuple{Vararg{Rule}})
    allkeys = Tuple(union(map(keys, rules)...))
    maxradii = Tuple(radius(rules, key) for key in allkeys)
    return NamedTuple{allkeys}(maxradii)
end
# Get radius of specific key from all rules
radius(rules::Tuple{Vararg{Rule}}, key::Symbol) =
    reduce(max, radius(ru) for ru in rules if key in keys(ru); init=0)

radius(rule::Rule, args...) = 0

"""
    ReturnRule <: Rule

Abstract supertype for rules that return value/s.

These must define methods of [`applyrule`](@ref).
"""
abstract type ReturnRule{R,W} <: Rule{R,W} end


"""
    Cellrule <: Rule

A `Rule` that only writes and uses a state from single cell of the read grids,
and has its return value written back to the same cell(s).

This limitation can be useful for performance optimisation,
such as wrapping rules in [`Chain`](@ref) so that no writes occur between rules.

`CellRule` is defined with :

```jldoctest yourcellrule
using DynamicGrids
struct MyCellRule{R,W} <: CellRule{R,W} end
# output
```

And applied as:

```jldoctest yourcellrule
function applyrule(data, rule::MyCellRule, state, I)
    state * 2
end
# output
applyrule (generic function with 1 method)
```

As the index `I` is provided in `applyrule`, you can use it to look up [`Aux`](@ref) data. 
"""
abstract type CellRule{R,W} <: ReturnRule{R,W} end

ruletype(::CellRule) = CellRule

"""
    Call <: CellRule

    Cell(f)
    Cell{R,W}(f)

A [`CellRule`](@ref) that applies a function `f` to the `R` grid value, 
or `Tuple` of values, and returns the `W` grid value or `Tuple` of values.

Especially convenient with `do` notation.

## Example

Double the cell value in grid `:a`:

```jldoctest Cell
using DynamicGrids
simplerule = Cell{:a}() do data, a, I
    2a
end
# output
Cell{:a,:a}(
    f = var"#1#2"
)
```

`data` is an [`AbstractSimData`](@ref) object, `a` is the cell value, and `I`
is a `Tuple` holding the cell index.

If you need to use multiple grids (a and b), use the `R`
and `W` type parameters. If you want to use external variables,
wrap the whole thing in a `let` block, for performance. This
rule sets the new value of `b` to the value of `a` to `b` times scalar `y`:

```jldoctest Cell
y = 0.7
rule = let y = y
    rule = Cell{Tuple{:a,:b},:b}() do data, (a, b), I
        a + b * y
    end
end
# output
Cell{Tuple{:a, :b},:b}(
    f = var"#3#4"{Float64}
)
```
"""
struct Cell{R,W,F} <: CellRule{R,W}
    "Function to apply to the R grid values"
    f::F
end
Cell{R,W}(; kw...) where {R,W} = _nofunctionerror(Cell)

@inline function applyrule(data, rule::Cell, read, I)
    let data=data, rule=rule, read=read, I=I
        rule.f(data, read, I)
    end
end

"""
    SetRule <: Rule

Abstract supertype for rules that manually write to the grid in some way.

These must define methods of [`applyrule!`](@ref).
"""
abstract type SetRule{R,W} <: Rule{R,W} end

ruletype(::SetRule) = SetRule

# Allow specifying CartesianIndex for `I`
applyrule!(data, rule, read, I) = applyrule(data, rule, read, CartesianIndex(I))

"""
    SetCellRule <: Rule

Abstract supertype for rules that can manually write to any cells of the
grid that they need to.

For example, `SetCellRule` is applied with like this, here simply adding 1 to the current cell:

```jldoctest SetCellRule
using DynamicGrids
struct MySetCellRule{R,W} <: SetCellRule{R,W} end

function applyrule!(data::AbstractSimData, rule::MySetCellRule{R,W}, state, I) where {R,W}
    # Add 1 to the cell 10 up and 10 accross
    I, isinbounds = inbounds(I .+ 10)
    isinbounds && add!(data[W], 1, I...)
    return nothing
end
# output
applyrule! (generic function with 1 method)
```

Note the `!` bang - this method alters the state of `data`.

To update the grid, you can use atomic operators [`add!`](@ref), [`sub!`](@ref),
[`min!`](@ref), [`max!`](@ref), and [`and!`](@ref), [`or!`](@ref) for `Bool`.
These methods safely combined writes from all grid cells - directly using `setindex!`
would cause bugs.

It there are multiple write grids, you will need to get the grid keys from
type parameters, here `W1` and `W2`:

```jldoctest SetCellRule
function applyrule(data, rule::MySetCellRule{R,Tuple{W1,W2}}, state, I) where {R,W1,W2}
    add!(data[W1], 1, I...)
    add!(data[W2], 2, I...)
    return nothing
end
# output
applyrule (generic function with 1 method)
```

DynamicGrids guarantees that:

- values written to anywhere on the grid do not affect other cells in
    the same rule at the same timestep.
- values written to anywhere on the grid are available to the next rule in the
    sequence, or in the next timestep if there are no remaining rules.
- if atomic operators like `add!` and `sub!` are always used to write to the grid,
    race conditions will not occur on any hardware.
"""
abstract type SetCellRule{R,W} <: SetRule{R,W} end

ruletype(::SetCellRule) = SetCellRule

"""
    SetCell <: SetCellRule

    SetCell(f)
    SetCell{R,W}(f)

A [`SetCellRule`](@ref) to manually write to the array where you need to.
`f` is passed a [`AbstractSimData`](@ref) object, the grid state or `Tuple` of grid 
states for the cell and a `Tuple{Int,Int}` index of the current cell.

To update the grid, you can use: [`add!`](@ref), [`sub!`](@ref) for `Number`,
and [`and!`](@ref), [`or!`](@ref) for `Bool`. These methods safely combined
writes from all grid cells - directly using `setindex!` would cause bugs.

## Example

Choose a destination cell and if it is in the grid, update it based on the 
state of both grids:

```jldoctest SetCell
using DynamicGrids
rule = SetCell{Tuple{:a,:b},:b}() do data, (a, b), I 
    dest = your_dest_pos_func(I)
    if isinbounds(data, dest)
        destval = your_dest_val_func(a, b)
        add!(data[:b], destval, dest...)
    end
end

# output
SetCell{Tuple{:a, :b},:b}(
    f = var"#1#2"
)
```
"""
struct SetCell{R,W,F} <: SetCellRule{R,W}
    "Function to apply to data, index and read grid values"
    f::F
end
SetCell{R,W}(; kw...) where {R,W} = _nofunctionerror(Set)

@inline function applyrule!(data, rule::SetCell, read, I)
    let data=data, rule=rule, read=read, I=I
        rule.f(data, read, I)
    end
end

"""
    NeighborhoodRule <: Rule

A Rule that only accesses a stencil centered around the current cell.
`NeighborhoodRule` is applied with the method:

```julia
applyrule(data::AbstractSimData, rule::MyNeighborhoodRule, state, I::Tuple{Int,Int})
```

`NeighborhoodRule` must have a `stencil` method or field, that holds
a [`Neighborhood`](@ref) object. `neighbors(rule)` returns an iterator
over the surrounding cell pattern defined by the `Neighborhood`.

For each cell in the grids the stencil buffer will be updated
for use in the `applyrule` method, managed to minimise array reads.

This allows memory optimisations and the use of high-perforance routines on the
stencil buffer. It also means that and no bounds checking is required in
stencil code.

For stencil rules with multiple read grids, the first is always
the one used for the stencil, the others are passed in as additional
state for the cell. Any grids can be written to, but only for the current cell.
"""
abstract type NeighborhoodRule{R,W} <: ReturnRule{R,W} end

ruletype(::NeighborhoodRule) = NeighborhoodRule

_val_ruletype(rule::T) where T<:Rule = Val{T}()

Stencils.neighbors(rule::NeighborhoodRule) = neighbors(stencil(rule))
Stencils.stencil(rule::NeighborhoodRule) = rule.stencil
Stencils.offsets(rule::NeighborhoodRule) = offsets(stencil(rule))
Stencils.cartesian_offsets(rule::NeighborhoodRule) = cartesian_offsets(stencil(rule))
Stencils.kernel(rule::NeighborhoodRule) = kernel(stencil(rule))
Stencils.kernelproduct(rule::NeighborhoodRule) = kernelproduct(stencil(rule))
Stencils.indices(rule::NeighborhoodRule, args...) = indices(stencil(rule), args...)
Stencils.distances(rule::NeighborhoodRule) = distances(stencil(rule))
Stencils.distance_zones(rule::NeighborhoodRule) = distance_zones(stencil(rule))
stencilkey(rule::NeighborhoodRule{R,W}) where {R,W} = R
# The first argument is for the stencil grid
stencilkey(rule::NeighborhoodRule{<:Tuple{R1,Vararg},W}) where {R1,W} = R1
Stencils.radius(rule::NeighborhoodRule, args...) = radius(stencil(rule))
Stencils.unsafe_neighbors(rule::NeighborhoodRule, A::Stencils.AbstractStencilArray, I::CartesianIndex) =
    Stencils.unsafe_neighbors(stencil(rule), A, I)
@inline function Stencils.rebuild(rule::NeighborhoodRule, neighbors::SVector)
    @set rule.stencil = Stencils.rebuild(stencil(rule), neighbors)
end
@inline function Stencils.unsafe_stencil(rule::NeighborhoodRule, A::AbstractArray, I::CartesianIndex)
    Stencils.rebuild(rule, unsafe_neighbors(rule, A, I))
end

"""
    Neighbors <: NeighborhoodRule

    Neighbors(f, stencil=Moor(1))
    Neighbors{R,W}(f, stencil=Moore())

A [`NeighborhoodRule`](@ref) that receives a [`Neighborhood`](@ref) object 
for the first `R` grid, followed by the cell value/s for the required grids, 
as with [`Cell`](@ref).

Returned value(s) are written to the `W` grid/s.

As with all [`NeighborhoodRule`](@ref), you do not have to check bounds at
grid edges, that is handled for you internally.

Using [`SparseOpt`](@ref) may improve stencil performance
when a specific value (often zero) is common and can be safely ignored.

## Example

Runs a game of life glider on grid `:a`:

```jldoctest
using DynamicGrids
const sum_states = (0, 0, 1, 0, 0, 0, 0, 0, 0), 
                   (0, 0, 1, 1, 0, 0, 0, 0, 0)
life = Neighbors{:a}(Moore(1)) do data, hood, a, I
    sum_states[a + 1][sum(hood) + 1]
end
init = Bool[
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 1 1 1 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
]
output = REPLOutput((; a=init); fps=25, tspan=1:50)
sim!(output, Life{:a}(); boundary=Wrap())
output[end][:a]

# output
5×15 Matrix{Bool}:
 0  0  1  0  1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  0  0  0  0  0  0  0  0  0  0

```
"""
struct Neighbors{R,W,F,N} <: NeighborhoodRule{R,W}
    "Function to apply to the stencil and R grid values"
    f::F
    "Defines the stencil of cells around the central cell"
    stencil::N
end
Neighbors{R,W}(; kw...) where {R,W} = _nofunctionerror(Neighbors)
Neighbors{R,W}(f; neighborhood=Moore(1), stencil=neighborhood) where {R,W} =
    Neighbors{R,W}(f, stencil)

@inline function applyrule(data, rule::Neighbors, read, I)
    let rule=rule, hood=stencil(rule), read=read, I=I
        rule.f(data, hood, read, I)
    end
end

"""
    Convolution <: NeighborhoodRule

    Convolution(kernel::AbstractArray)
    Convolution{R,W}(kernel::AbstractArray)

A [`NeighborhoodRule`](@ref) that runs a convolution kernel over the grid.

`kernel` must be a square matrix.

## Performance

Small radius convolutions in DynamicGrids.jl will be comparable or even faster than using
DSP.jl or ImageConvolutions.jl. As the radius increases these packages will be a lot faster.

But `Convolution` is convenient to chain into a simulation, and combined with some other
rules. It should perform reasonably well with all but very large kernels.

Values are not normalised, so make sure the kernel sums to `1` if you need that.

## Example

A streaking convolution that looks a bit like sand blowing.

Swap out the matrix values to change the pattern.

```julia 
using DynamicGrids, DynamicGridsGtk
streak = Convolution([0.0 0.01 0.48; 
                      0.0 0.5  0.01; 
                      0.0 0.0  0.0])
output = GtkOutput(rand(500, 500); tspan = 1:1000, fps=100)
sim!(output, streak; boundary=Wrap())
```
"""
struct Convolution{R,W,N} <: NeighborhoodRule{R,W}
    "The stencil of cells around the central cell"
    stencil::N
    Convolution{R,W}(stencil::N) where {R,W,N<:Stencil} = new{R,W,N}(stencil)
end
Convolution{R,W}(A::AbstractArray) where {R,W} = Convolution{R,W}(Kernel(SArray{Tuple{size(A)...}}(A)))
Convolution{R,W}(; stencil) where {R,W} = Convolution{R,W}(stencil)

@inline applyrule(data, rule::Convolution, read, I) = kernelproduct(stencil(rule))

"""
    SetNeighborhoodRule <: SetRule

A [`SetRule`](@ref) that only writes to its stencil, and does not need to bounds-check.

[`indices`](@ref) and [`offsets`](@ref) are useful iterators for modifying
stencil values. 

`SetNeighborhoodRule` rules must return a [`Neighborhood`](@ref) object from the function 
`stencil(rule)`. By default this is `rule.stencil`. If this property exists, 
no interface methods are required.
"""
abstract type SetNeighborhoodRule{R,W} <: SetRule{R,W} end

ruletype(::SetNeighborhoodRule) = SetNeighborhoodRule

Stencils.stencil(rule::SetNeighborhoodRule) = rule.stencil
Stencils.offsets(rule::SetNeighborhoodRule) = offsets(stencil(rule))
Stencils.kernel(rule::SetNeighborhoodRule) = kernel(stencil(rule))
Stencils.indices(rule::SetNeighborhoodRule, args...) = indices(stencil(rule), args...)
Stencils.radius(rule::SetNeighborhoodRule, args...) = radius(stencil(rule))
stencilkey(rule::SetNeighborhoodRule{R,W}) where {R,W} = R
stencilkey(rule::SetNeighborhoodRule{<:Tuple{R1,Vararg},W}) where {R1,W} = R1

"""
    SetNeighbors <: SetNeighborhoodRule

    SetNeighbors(f, stencil=Moor(1))
    SetNeighbors{R,W}(f, stencil=Moor(1))

A [`SetCellRule`](@ref) to manually write to the array with the specified
stencil. Indexing outside the stencil is undefined behaviour.

Function `f` is passed four arguments: a [`SimData`](@ref) object, the specified
[`Neighborhood`](@ref) object, the grid state or `Tuple` of grid states for the cell, 
and the `Tuple{Int,Int}` index of the current cell.

To update the grid, you can use: [`add!`](@ref), [`sub!`](@ref) for `Number`,
and [`and!`](@ref), [`or!`](@ref) for `Bool`. These methods can be safely combined
writes from all grid cells.

Directly using `setindex!` is possible, but may cause bugs as multiple cells
may write to the same location in an unpredicatble order. As a rule, directly
setting a stencil index should only be done if it always sets the samevalue -
then it can be guaranteed that any writes from othe grid cells reach the same result.

[`neighbors`](@ref), [`offsets`](@ref) and [`indices`](@ref) are useful methods for
`SetNeighbors` rules.

## Example

This example adds a value to all neighbors:

```jldoctest
using DynamicGrids

rule = SetNeighbors{:a}() do data, stencil, a, I
    add_to_neighbors = your_func(a)
    for pos in indices(stencil)
        add!(data[:b], add_to_neighbors, pos...)
    end
end
# output
SetNeighbors{:a,:a}(
    f = var"#1#2"
    stencil = Moore{1, 2, 8, Nothing}
)
```
"""
struct SetNeighbors{R,W,F,N} <: SetNeighborhoodRule{R,W}
    "Function to apply to the data, index and R grid values"
    f::F
    "The stencil of cells around the central cell"
    stencil::N
end
SetNeighbors{R,W}(; kw...) where {R,W} = _nofunctionerror(SetNeighbors)
SetNeighbors{R,W}(f; neighborhood=Moore(1), stencil=neighborhood) where {R,W} = 
    SetNeighbors{R,W}(f, stencil)

@inline function applyrule!(data, rule::SetNeighbors, read, I)
    let data=data, hood=stencil(rule), rule=rule, read=read, I=I
        rule.f(data, hood, read, I)
    end
end

"""
    SetGridRule <: Rule

A `Rule` applies to whole grids. This is used for operations that don't benefit from
having stencil buffering or looping over the grid handled for them, or any specific
optimisations. Best suited to simple functions like `rand!(grid)` or using convolutions
from other packages like DSP.jl. They may also be useful for doing other custom things that
don't fit into the DynamicGrids.jl framework during the simulation.

Grid rules specify the grids they want and are sequenced just like any other grid.

```julia
struct MySetGridRule{R,W} <: SetGridRule{R,W} end
```

And applied as:

```julia
function applyrule!(data::AbstractSimData, rule::MySetGridRule{R,W}) where {R,W}
    rand!(data[W])
end
```
"""
abstract type SetGridRule{R,W} <: Rule{R,W} end

ruletype(::SetGridRule) = SetGridRule

"""
    SetGrid{R,W}(f)

Apply a function `f` to fill whole grid/s.

Broadcasting is a good way to update values.

## Example

This example simply sets grid `a` to equal grid `b`:

```jldoctest
using DynamicGrids
rule = SetGrid{:a,:b}() do a, b
    b .= a
end

# output
SetGrid{:a,:b}(
    f = var"#1#2"
)
```
"""
struct SetGrid{R,W,F} <: SetGridRule{R,W}
    "Function to apply to the read values"
    f::F
end
SetGrid{R,W}(; kw...) where {R,W} = _nofunctionerror(SetGrid)

@inline function applyrule!(data, rule::SetGrid{R,W}) where {R,W}
    let data = data, rule=rule
        read = if R isa Symbol
            gridview(getindex(data, R))
        else
            map(w -> gridview(getindex(data, w)), _asiterable(R))
        end
        write = if W isa Symbol
            gridview(getindex(data, W))
        else
            map(w -> gridview(getindex(data, w)), _asiterable(W))
        end
        rule.f(data, read, write)
    end
end

# Utils

@noinline _nofunctionerror(t) =
    throw(ArgumentError("No function passed to $t. Did you mean to use a do block?"))

# Check number of args passed in as we don't get a normal method
# error because of the splatted args in the default constructor.
# @generated function _checkfields(::Type{T}, args::A) where {T,A<:Tuple}
    # length(fieldnames(T)) == length(fieldnames(A)) ? :(nothing) : :(_fielderror(T, args))
# end

# @noinline function _fielderror(T, args)
    # throw(ArgumentError("$T has $(length(fieldnames(T))) fields: $(fieldnames(T)), you have used $(length(args))"))
# end


