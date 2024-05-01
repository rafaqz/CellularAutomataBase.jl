"""
    ParameterSource

Abstract supertypes for parameter source wrappers such as [`Aux`](@ref),
[`Grid`](@ref) and [`Delay`](@ref). These allow flexibility in that parameters
can be retreived from various data sources not specified when the rule is written.
"""
abstract type ParameterSource end

"""
    get(data::AbstractSimData, source::ParameterSource, I...)
    get(data::AbstractSimData, source::ParameterSource, I::Union{Tuple,CartesianIndex})

Allows parameters to be taken from a single value or a [`ParameterSource`](@ref)
such as another [`Grid`](@ref), an [`Aux`](@ref) array, or a [`Delay`](@ref).

Other `source` objects are used as-is without indexing with `I`.
"""
@propagate_inbounds Base.get(data::AbstractSimData, val, I...) = val
@propagate_inbounds Base.get(data::AbstractSimData, key::ParameterSource, i::Integer, I::Integer...) =
    get(data, key, (i, I...))
@propagate_inbounds Base.get(data::AbstractSimData, key::ParameterSource, I::CartesianIndex) =
    get(data, key, Tuple(I))

"""
    Aux <: ParameterSource

    Aux{K}()
    Aux(K::Symbol)

Use auxilary array with key `K` as a parameter source.

Implemented in rules with:

```julia
get(data, rule.myparam, I)
```

When an `Aux` param is specified at rule construction with:

```julia
rule = SomeRule(; myparam=Aux{:myaux})
output = ArrayOutput(init; aux=(myaux=myauxarray,))
```

If the array is a DimensionalData.jl `DimArray` with a `Ti` (time)
dimension, the correct interval will be selected automatically,
precalculated for each timestep so it has no significant overhead.

Currently this is cycled by default. Note that cycling may be incorrect
when the simulation timestep (e.g. `Week`) does not fit
equally into the length of the time dimension (e.g. `Year`).
This will reuire a `Cyclic` index mode in DimensionalData.jl in future
to correct this problem.
"""
struct Aux{K} <: ParameterSource end
Aux(key::Symbol) = Aux{key}()
Aux(key::Val{K}) where K = Aux{K}()

_unwrap(::Aux{X}) where X = X
_unwrap(::Type{<:Aux{X}}) where X = X

@inline aux(nt::NamedTuple, ::Aux{Key}) where Key = nt[Key]

@propagate_inbounds function Base.get(data::AbstractSimData, key::Aux, I::Tuple)
    _getaux(data, key, I)
end
@propagate_inbounds function Base.get(data::AbstractSimData, key::Aux)
    aux(data, key)
end

# _getaux
# Get the value of an auxilliary array at index I and/or the synchonised time-step
@propagate_inbounds function _getaux(
    data::AbstractSimData, key::Union{Aux,Symbol}, I::Tuple
)
    _getaux(aux(data, key), data, key, I)
end
# For an Array just return the value for the index
@propagate_inbounds function _getaux(
    A::AbstractArray, data::AbstractSimData, key::Union{Aux,Symbol}, I::Tuple
)
    A[I...]
end
# For a DimArray with a time dimension we return the value at the
# current auxframe, also using the index `I` if aux is multidimensional.
@propagate_inbounds function _getaux(
    A::AbstractDimArray{<:Any,1}, data::AbstractSimData, key::Union{Aux,Symbol}, I::NTuple{N}
) where N
    hasdim(A, TimeDim) ? A[auxframe(data, key)] : A[I...]
end
@propagate_inbounds function _getaux(
    A::AbstractDimArray{<:Any,N1}, data::AbstractSimData, key::Union{Aux,Symbol}, I::NTuple{N2}
) where {N1,N2}
    if hasdim(A, TimeDim)
        last(dims(A)) isa TimeDim || throw(ArgumentError("The time dimensions in aux data must be the last dimension"))
        A[ntuple(i -> I[i], Val{N1-1}())..., auxframe(data, key)]
    else
        A[ntuple(i -> I[i], Val{N1-1}())...]
    end
end


# boundscheck_aux
# Bounds check the aux arrays ahead of time
function boundscheck_aux(data::AbstractSimData, key::Aux)
    boundscheck_aux(data, aux(data), key)
end
function boundscheck_aux(data::AbstractSimData, A::Nothing, key::Aux{Key}) where Key
    _auxmissingerror(Key, a)
end
function boundscheck_aux(data::AbstractSimData, aux::NamedTuple, key::Aux{Key}) where Key
    boundscheck_aux(data, aux[_unwrap(key)], key)
end
function boundscheck_aux(data::AbstractSimData, A::AbstractArray, key::Aux{Key}) where Key
    size(A) === size(data) || _auxsizeerror(Key, size(A), size(data))
end
function boundscheck_aux(data::AbstractSimData, A::AbstractDimArray{<:Any,1}, key::Aux{Key}) where Key
    hasdim(A, TimeDim) || size(data) === size(A) || _auxsizeerror(Key, size(A), size(data))
end
function boundscheck_aux(data::AbstractSimData, A::AbstractDimArray, key::Aux{Key}) where Key
    if hasdim(A, TimeDim)
        asize = size(otherdims(A, TimeDim))
        asize == size(data) || _auxsizeerror(Key, asize, size(data))
    else
        size(A) == size(data) || _auxsizeerror(Key, size(A), size(data))
    end
end

# _calc_auxframe
# Calculate the frame to use in the aux data for this timestep.
# This uses the index of any AbstractDimArray, which must be a
# matching type to the simulation tspan.
# This is called from _updatetime in simulationdata.jl
_calc_auxframe(data::AbstractSimData) = _calc_auxframe(aux(data), data)
function _calc_auxframe(aux::NamedTuple{K}, data::AbstractSimData) where K
    map((A, k) -> _calc_auxframe(A, data, k), aux, NamedTuple{K}(K))
end
function _calc_auxframe(A::AbstractDimArray, data, key)
    hasdim(A, TimeDim) || return nothing
    timedim = dims(A, TimeDim)
    curtime = currenttime(data)
    if !hasselection(timedim, Contains(curtime)) 
        if lookup(timedim) isa Cyclic
            if sampling(timedim) isa Points
                throw(ArgumentError("$(_no_valid_time(timedim,key, curtime)) Did you mean to use `Intervals` for the time dimension `sampling`? `Contains` on `Points` defaults to `At`, and must be exact."))
            else
                throw(ArgumentError("$(_no_valid_time(timedim,key, curtime))"))
            end
        elseif sampling(timedim) isa Points
            throw(ArgumentError("$(_no_valid_time(timedim,key, curtime)) Did you mean to use `Intervals` for the time dimension `sampling`? `Contains` on `Points` defaults to `At`, and must be exact."))
        else
            throw(ArgumentError("$(_no_valid_time(timedim,key, curtime)) Did you mean to use a `Cyclic` lookup for the time dimension of the array?"))
        end
    end
    return DimensionalData.selectindices(timedim, Contains(curtime))
end
_calc_auxframe(args...) = nothing

_no_valid_time(timedim, key, curtime) = "Time dimension over $(bounds(timedim)) of aux `$key` has no valid selection for `Contains($curtime)`."

"""
    Grid <: ParameterSource

    Grid{K}()
    Grid(K::Symbol)

Use grid with key `K` as a parameter source.

Implemented in rules with:

```julia
get(data, rule.myparam, I)
```

And specified at rule construction with:

```julia
SomeRule(; myparam=Grid(:somegrid))
```
"""
struct Grid{K} <: ParameterSource end
Grid(key::Symbol) = Grid{key}()

_unwrap(::Grid{X}) where X = X
_unwrap(::Type{<:Grid{X}}) where X = X

@propagate_inbounds Base.get(data::AbstractSimData, key::Grid{K}, I::Tuple) where K =
    data[K][I...]


"""
    AbstractDelay <: ParameterSource

Abstract supertype for [`ParameterSource`](@ref)s that use data from a grid
with a time delay.

$EXPERIMENTAL
"""
abstract type AbstractDelay{K} <: ParameterSource end

@inline frame(delay::AbstractDelay) = delay.frame

@propagate_inbounds function Base.get(
    data::AbstractSimData, delay::AbstractDelay{K}, I::Tuple
) where K
    _getdelay(frames(data)[frame(delay, data)], delay, I)
end

# The output may store a single Array or a NamedTuplea.
@propagate_inbounds function _getdelay(frame::AbstractArray, ::AbstractDelay, I::Tuple)
    frame[I...]
end
@propagate_inbounds function _getdelay(
    frame::NamedTuple, delay::AbstractDelay{K}, I::Tuple
) where K
    _getdelay(frame[K], delay, I)
end

"""
    Delay <: AbstractDelay

    Delay{K}(steps)

`Delay` allows using a [`Grid`](@ref) from previous timesteps as a parameter source as
a field in any `Rule` that uses `get` to retrieve it's parameters.

It must be coupled with an output that stores all frames, so that `@assert
DynamicGrids.isstored(output) == true`.  With [`GraphicOutput`](@ref)s this may be
acheived by using the keyword argument `store=true` when constructing the output object.

# Type Parameters

- `K::Symbol`: matching the name of a grid in `init`.

# Arguments

- `steps`: As a user supplied parameter, this is a multiple of the step size of the output
    `tspan`. This is automatically replaced with an integer for each step. Used within the
    code in a rule, it must be an `Int` number of frames, for performance.

# Example

```julia
SomeRule(;
    someparam=Delay(:grid_a, Month(3))
    otherparam=1.075
)
```

$EXPERIMENTAL
"""
struct Delay{K,S} <: AbstractDelay{K}
    steps::S
end
Delay{K}(steps::S) where {K,S} = Delay{K,S}(steps)

ConstructionBase.constructorof(::Delay{K}) where K = Delay{K}

steps(delay::Delay) = delay.steps

# _to_frame
# Replace the delay step size with an integer for fast indexing
# checking that the delay matches the simulation step size.
# Delays at the start just use the init frame.
function _to_frame(delay::Delay{K}, data::AbstractSimData) where K
    nsteps = steps(delay) / step(tspan(data))
    isteps = trunc(Int, nsteps)
    nsteps == isteps || _delaysteperror(delay, step(tspan(data)))
    frame = max(currentframe(data) - isteps, 1)
    Frame{K}(frame)
end

"""
    Lag <: AbstractDelay

    Lag{K}(frames::Int)

`Lag` allows using a [`Grid`](@ref) from a specific previous frame from within a rule,
using `get`. It is similar to [`Delay`](@ref), but an integer amount of frames should be
used, instead of a quantity related to the simulation `tspan`. The lower bound is the first
frame.

# Type Parameter

- `K::Symbol`: matching the name of a grid in `init`.

# Argument

- `frames::Int`: number of frames to lag by, 1 or larger.

# Example

```julia
SomeRule(;
    someparam=Lag(:grid_a, Month(3))
    otherparam=1.075
)
```

$EXPERIMENTAL
"""
struct Lag{K} <: AbstractDelay{K}
    nframes::Int
end

# convert the Lag to a Frame object
# Lags at the start just use the init frame.
function _to_frame(ps::Lag{K}, data::AbstractSimData) where K
    Frame{K}(frame(ps, data))
end

@inline frame(ps::Lag, data) = max(1, currentframe(data) - ps.nframes)

"""
    Frame <: AbstractDelay

    Frame{K}(frame)

`Frame` allows using a [`Grid`](@ref) from a specific previous timestep from within
a rule, using `get`. It should only be used within rule code, not as a parameter.

# Type Parameter

- `K::Symbol`: matching the name of a grid in `init`.

# Argument

- `frame::Int`: the exact frame number to use.

$EXPERIMENTAL
"""
struct Frame{K} <: AbstractDelay{K}
    frame::Int
    Frame{K}(frame::Int) where K = new{K}(frame)
end

ConstructionBase.constructorof(::Frame{K}) where K = Frame{K}

@inline frame(ps::Frame) = ps.frame
@inline frame(ps::Frame, data) = frame(ps)


# Delay utils

# Types to ignore when flattening rules to <: AbstractDelay
const DELAY_IGNORE = Union{Function,SArray,AbstractDict,Number}

@inline function hasdelay(rules::Tuple)
    # Check for Delay as parameter or used in rule code
    length(_getdelays(rules)) > 0 || any(map(needsdelay, rules))
end

# needsdelay
# Speficy that a rule needs a delay frames present
# to run. This will throw an early error if the Output
# does not store frames, instead of an indexing error during
# the simulation.
needsdelay(rule::Rule) = false

# _getdelays
# Get all the delays found in fields of the rules tuple
_getdelays(rules::Tuple) = Flatten.flatten(rules, AbstractDelay, DELAY_IGNORE)

# _setdelays
# Update any AbstractDelay anywhere in the rules Tuple.
# These are converted to Frame objects so the calculation
# happens only once for each timestep, instead of for each cell.
function _setdelays(rules::Tuple, data::AbstractSimData)
    delays = _getdelays(rules)
    if length(delays) > 0
        newdelays = map(d -> _to_frame(d, data), delays)
        _setdelays(rules, newdelays)
    else
        rules
    end
end
_setdelays(rules::Tuple, delays::Tuple) =
    Flatten.reconstruct(rules, delays, AbstractDelay, DELAY_IGNORE)

# Errors
@noinline _notstorederror() =
    throw(ArgumentError("Output does not store frames, which is needed for a Delay. Use ArrayOutput or any GraphicOutput with `store=true` keyword"))
@noinline _delaysteperror(delay::Delay{K}, simstep) where K =
    throw(ArgumentError("Delay $K size $(steps(delay)) is not a multiple of simulations step $simstep"))
@noinline _auxmissingerror(key, aux) = error("Aux data $key is not present in aux")
@noinline _auxsizeerror(key, asize, gsize) =
    error("Aux data $key is size $asize does not match grid size $gsize")
