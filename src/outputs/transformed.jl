"""
    TransformedOutput(f, init; tspan::AbstractRange, kw...) 

An output that stores the result of some function `f` of the grid/s.

# Arguments

- `f`: a function or functor that accepts an `AbstractArray` or `NamedTuple` of
    `AbstractArray` with names matching `init`. The `AbstractArray` will be a view into 
    the grid the same size as the init grids, removing any padding that has been added.
- `init`: initialisation `Array` or `NamedTuple` of `Array`

# Keywords

- `tspan`: `AbstractRange` timespan for the simulation
- `aux`: NamedTuple of arbitrary input data. Use `get(data, Aux(:key), I...)` 
    to access from a `Rule` in a type-stable way.
- `mask`: `BitArray` for defining cells that will/will not be run.
- `padval`: padding value for grids with stencil rules. The default is `zero(eltype(init))`.
- `eltype`: force the output vector to have a specific element type.
    Useful if `f` doesn't always return the same type.

$EXPERIMENTAL
"""
mutable struct TransformedOutput{T,F,A<:AbstractVector{T},E,B} <: Output{T,A} 
    f::F
    frames::A
    running::Bool
    extent::E
    buffer::B
end
function TransformedOutput(f::Function, init::Union{NamedTuple,AbstractMatrix}; 
    extent=nothing, 
    eltype=nothing,
    tspan, 
    kw...
)
    _check_grids_are_isbits(init)
    # We have to handle some things manually as we are changing the standard output frames
    extent = extent isa Nothing ? Extent(; init=init, tspan, kw...) : extent
    # Define buffers to copy to before applying `f`
    buffer = _replicate_init(init, replicates(extent))
    f1 = f(buffer)
    if buffer isa NamedTuple
        map(buffer) do b
            b .= (zero(first(b)),)
        end
    else
        buffer .= (zero(first(buffer)),)
    end
    # Build simulation frames from the output of `f` for empty frames
    frames = if isnothing(eltype)
        [deepcopy(f1) for _ in eachindex(tspan)]
    else
        eltype[deepcopy(f1) for _ in eachindex(tspan)]
    end
    # Set the first frame to the output of `f` for `init`
    frames[1] = f1

    return TransformedOutput(f, frames, false, extent, buffer)
end
function TransformedOutput(init; kw...)
    throw(ArgumentError("TransformedOutput must be passed a function and the init grid(s) as arguments"))
end

function storeframe!(o::TransformedOutput, data::AbstractSimData) 
    transformed = _transform_grids(o, grids(data))
    i = frameindex(o, data) 
    # Copy the transformed grid/s to the output frames, 
    # instead of just assigning (see issue #169)
    o[i] = _copytransformed!(o[i], transformed)
end


# Multi/named grid simulation, f is passed a NamedTuple
_transform_grids(o::TransformedOutput, grids::NamedTuple) = o.f(grids)
# Single unnamed grid simulation, f is passed an AbstractArray
_transform_grids(o::TransformedOutput, grids::NamedTuple{(DEFAULT_KEY,)}) = o.f(first(grids))

# Copy arrays manually as reducing functions can return the original object without copy.
_copytransformed!(dest::NamedTuple, src::NamedTuple) = map(_copytransformed!, dest, src)
_copytransformed!(dest::AbstractArray, src::AbstractArray) = dest .= src
# Non-array output is just assigned
_copytransformed!(dest, src) = src
_copytransformed!(dest::StaticArray, src::StaticArray) = src

init_output_grids!(o::TransformedOutput, init) = nothing

Adapt.adapt_structure(to, o::TransformedOutput) = o
