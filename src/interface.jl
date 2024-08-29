"""
    applyrule(data::AbstractSimData, rule::Rule, state, index::NTuple{N,Int}) -> new cell value/s
    applyrule(data::AbstractSimData, rule::Rule, state, index::CaresianIndex) -> new cell value/s
    applyrule(rule::Rule, state, index::NTuple{N,Int}) -> new cell value/s
    applyrule(rule::Rule, state, index::CaresianIndex) -> new cell value/s
    applyrule(rule::Rule, state) -> new cell value/s
    applyrule(rule::Rule) -> new cell value/s

Apply a rule to the cell state and return values to write to the grid(s).

This is called in `maprule!` methods during the simulation,
not by the user. Custom `Rule` implementations must define this method,
but may use any of the combinations of arguments show above.

## Arguments

- `data` : [`AbstractSimData`](@ref)
- `rule` : [`Rule`](@ref)
- `state`: the value(s) of the current cell
- `index`: a (row, column) tuple of `Int` for the current cell coordinates

Returns the value/s to be written to the current cell of
the grids specified by the `W` type parameter. 

For a single grid, return a single value. For multiple grids, a `Tuple`
with values for each grid in the order specified in `W`.
"""
function applyrule end

"""
    applyrule!(data::AbstractSimData, rule::{R,W}, state, index::Tuple{Int,Int}) -> Nothing

Apply a rule to the cell state and manually write to the grid data array.
Used in all rules inheriting from [`SetCellRule`](@ref).

This is called in internal `maprule!` methods during the simulation, not by
the user. Custom [`SetCellRule`](@ref) implementations must define this method.

Only grids specified with the `W` type parameter will be writable from `data`.

## Arguments

- `data` : [`AbstractSimData`](@ref)
- `rule` : [`Rule`](@ref)
- `state`: the value(s) of the current cell
- `index`: a (row, column) tuple of Int for the current cell coordinates - `t`: the current time step
"""
function applyrule! end

"""
    modifyrule(rule::Rule, data::AbstractSimData) -> Rule

Precalculates rule fields at each timestep. Define this method if a [`Rule`](@ref)
has fields that need to be updated over time.

`Rule`s are immutable (it's faster and works on GPU), so `modifyrule` is
expected to return a new rule object with changes applied to it. Setfield.jl or
Acessors.jl may help with updating the immutable struct.

The default behaviour is to return the existing rule without change. Updated rules
are discarded after use, and the `rule` argument is always the original object passed in.

# Example

We define a rule with a parameter that is the total sum of the grids current,
and update it for each time-step using `modifyrule`.

This could be used to simulate top-down control e.g. a market mechanism in a
geographic model that includes agricultural economics.

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
"""
function modifyrule end

"""
    add!(data::WritableGridData, x, I...)

Add the value `x` to a grid cell.

## Example useage

```jldoctest
using DynamicGrids
rule = SetCell{:a}() do data, a, cellindex
    dest, is_inbounds = inbounds(data, (jump .+ cellindex)...)

    # Update spotted cell if it's on the grid
    is_inbounds && add!(data[:a], state, dest...)
end

# output
SetCell{:a,:a}(
    f = var"#1#2"
)
```
"""
function add! end

"""
    sub!(data::WritableGridData, x, I...)

Subtract the value `x` from a grid cell. See `add!` for example usage.
"""
function sub! end

"""
    min!(data::WritableGridData, x, I...)

Set a gride cell to the minimum of `x` and the current value. See `add!` for example usage.
"""
function min! end

"""
    max!(data::WritableGridData, x, I...)

Set a gride cell to the maximum of `x` and the current value. See `add!` for example usage.
"""
function max! end

"""
    and!(data::WritableGridData, x, I...)
    and!(A::AbstractArray, x, I...)

Set the grid cell `c` to `c & x`. See `add!` for example usage.
"""
function and! end

"""
    or!(data::WritableGridData, x, I...)
    or!(A::AbstractArray, x, I...)

Set the grid cell `c` to `c | x`. See `add!` for example usage.
"""
function or! end

"""
    xor!(data::WritableGridData, x, I...)
    xor!(A::AbstractArray, x, I...)

Set the grid cell `c` to `xor(c, x)`. See `add!` for example usage.
"""
function xor! end

"""
    inbounds(data::AbstractSimData, I::Tuple) -> Tuple{NTuple{2,Int}, Bool}
    inbounds(data::AbstractSimData, I...) -> Tuple{NTuple{2,Int}, Bool}

Check grid boundaries for a coordinate before writing in [`SetCellRule`](@ref).

Returns a `Tuple` containing a coordinates `Tuple` and a `Bool` - `true`
if the cell is inside the grid bounds, `false` if not.

[`BoundaryCondition`](@ref) of type [`Remove`](@ref) returns the coordinate and `false` 
to skip coordinates that boundary outside of the grid.

[`Wrap`](@ref) returns a tuple with the current position or it's
wrapped equivalent, and `true` as it is allways in-bounds.

[`Reflect`](@ref) returns a tuple with the current position or it's
reflected equivalent, and `true` as it is allways in-bounds.
"""
function inbounds end

"""
    isinbounds(data, I::Tuple) -> Bool
    isinbounds(data, I...) -> Bool

Check that a coordinate is within the grid, usually in [`SetCellRule`](@ref).

Unlike [`inbounds`](@ref), [`BoundaryCondition`](@ref) status is ignored.
"""
function isinbounds end

"""
    init(obj) -> Union{AbstractArray,NamedTUple}

Retrieve the mask from an [`Output`](@ref), [`Extent`](@ref) or [`AbstractSimData`](@ref) object.
"""
function init end

"""
    mask(obj) -> AbstractArray

Retrieve the mask from an [`Output`](@ref), [`Extent`](@ref) or [`AbstractSimData`](@ref) object.
"""
function mask end

"""
    aux(obj, [key])

Retrieve auxilary data `NamedTuple` from an [`Output`](@ref),
[`Extent`](@ref) or [`AbstractSimData`](@ref) object.

Given `key` specific data will be returned. `key` should be a
`Val{:symbol}` for type stability and zero-cost access inside rules.
`Symbol` will also work, but may be slow.
"""
function aux end

"""
    tspan(obj) -> AbstractRange

Retrieve the time-span `AbstractRange` from an [`Output`](@ref),
[`Extent`](@ref) or [`AbstractSimData`](@ref) object.
"""
function tspan end

"""
    timestep(obj)

Retrieve the timestep size from an [`Output`](@ref),
[`Extent`](@ref), [`Ruleset`](@ref) or [`AbstractSimData`](@ref) object.

This will be in whatever type/units you specify in `tspan`.
"""
function timestep end

"""
    currentframe(simdata::AbstractSimData) -> Int

Retrieve the current simulation frame a [`AbstractSimData`](@ref) object.
"""
function currentframe end

"""
    currenttime(simdata::AbstractSimData)

Retrieve the current simulation time from a [`AbstractSimData`](@ref) object.

This will be in whatever type/units you specify in `tspan`.
"""
function currenttime end

"""
    currenttimestep(simdata::AbstractSimData)

Retrieve the current timestep from a [`AbstractSimData`](@ref) object.

This may be different from the `timestep`. If the timestep is `Month`,
`currenttimestep` will return `Seconds` for the length of the specific month.
"""
function currenttimestep end

