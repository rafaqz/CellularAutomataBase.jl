
const GRAPHICCONFIG_KEYWORDS = """
- `fps::Real`: Frames per second.
- `store::Bool`: Whether to store frames like `ArrayOutput` or to disgard
    them after visualising. Very long simulation runs may fill available 
    memory when `store=true`.
"""

"""
    GraphicConfig

    GraphicConfig(; fps=25.0, store=false)

Config and variables for graphic outputs.

# Keywords

$GRAPHICCONFIG_KEYWORDS
"""
mutable struct GraphicConfig{FPS,TS}
    fps::FPS
    store::Bool
    timestamp::TS
    stampframe::Int
    stoppedframe::Int
end
GraphicConfig(; fps=25.0, store=false, kw...) = GraphicConfig(fps, store, 0.0, 1, 1)

fps(gc::GraphicConfig) = gc.fps
store(gc::GraphicConfig) = gc.store
timestamp(gc::GraphicConfig) = gc.timestamp
stampframe(gc::GraphicConfig) = gc.stampframe
stoppedframe(gc::GraphicConfig) = gc.stoppedframe
setfps!(gc::GraphicConfig, x) = gc.fps = x
function settimestamp!(o::GraphicConfig, frame::Int)
    o.timestamp = time()
    o.stampframe = frame
end

setstoppedframe!(gc::GraphicConfig, frame::Int) = gc.stoppedframe = frame

const GRAPHICOUTPUT_KEYWORDS = """
## [`Extent`](@ref) keywords:
$EXTENT_KEYWORDS

An `Extent` object can be also passed to the `extent` keyword, and other keywords will be ignored.

## [`GraphicConfig`](@ref) keywords:
$GRAPHICCONFIG_KEYWORDS

A `GraphicConfig` object can be also passed to the `graphicconfig` keyword, and other keywords will be ignored.
"""

"""
    GraphicOutput <: Output

Abstract supertype for [`Output`](@ref)s that display the simulation frames.

All `GraphicOutputs` must have a [`GraphicConfig`](@ref) object
and define a [`showframe`](@ref) method.

See [`REPLOutput`](@ref) for an example.

# User Arguments for all `GraphicOutput`:

- `init`: initialisation `AbstractArray` or `NamedTuple` of `AbstractArray`

# Minimum user keywords for all `GraphicOutput`:

$GRAPHICOUTPUT_KEYWORDS

## Internal keywords for constructors of objects extending `GraphicOutput`: 

The default constructor will generate these objects and pass them to the inheriting 
object constructor, which must accept the following keywords:

- `frames`: a `Vector` of simulation frames (`NamedTuple` or `Array`). 
- `running`: A `Bool`.
- `extent` an [`Extent`](@ref) object.
- `graphicconfig` a [`GraphicConfig`](@ref)object.

Users can also pass in these entire objects if required.
"""
abstract type GraphicOutput{T,F} <: Output{T,F} end

# Generic GraphicOutput constructor. Converts an init array to vector of arrays.
function (::Type{T})(
    init::Union{NamedTuple,AbstractArray}; 
    extent=nothing, graphicconfig=nothing, store=nothing, kw...
) where T <: GraphicOutput
    extent = extent isa Nothing ? Extent(; init=init, kw...) : extent
    store = check_stored(extent, store)
    graphicconfig = graphicconfig isa Nothing ? GraphicConfig(; store, kw...) : graphicconfig
    frames = [_replicate_init(init, replicates(extent))]
    T(; frames, running=false, extent, graphicconfig, store, kw...)
end

function check_stored(extent, store) 
    if ndims(extent) == 1 && !(store === true)
        if store === false
            @warn "Setting `store=false` may cause errors when visualising a 1-dimensional simulation"
        else
            store = true
        end
    else
        store = store isa Nothing ? false : store
    end
    return store
end
 
graphicconfig(o::Output) = GraphicConfig()
graphicconfig(o::GraphicOutput) = o.graphicconfig

# Forward getters and setters to GraphicConfig object ####################################
fps(o::GraphicOutput) = fps(graphicconfig(o))
timestamp(o::GraphicOutput) = timestamp(graphicconfig(o))
stampframe(o::GraphicOutput) = stampframe(graphicconfig(o))
stoppedframe(o::GraphicOutput) = stoppedframe(graphicconfig(o))
isstored(o::GraphicOutput) = store(o)
store(o::GraphicOutput) = store(graphicconfig(o))

setfps!(o::GraphicOutput, x) = setfps!(graphicconfig(o), x)
settimestamp!(o::GraphicOutput, f) = settimestamp!(graphicconfig(o), f)
setstoppedframe!(o::GraphicOutput, f) = setstoppedframe!(graphicconfig(o), f)

# Delay output to maintain the frame rate
maybesleep(o::GraphicOutput, f) =
    sleep(max(0.0, timestamp(o) + (f - stampframe(o))/fps(o) - time()))
isshowable(o::GraphicOutput, f) = true

# Store frames, and show them graphically ################################################
function storeframe!(o::GraphicOutput, data)
    f = frameindex(o, data)
    if f > length(o)
        _pushgrid!(eltype(o), o)
    end
    if isstored(o)
        _storeframe!(o, eltype(o), data)
    end
    if isshowable(o, currentframe(data))
        showframe(o, data)
    end
    return nothing
end

# Add additional grids if they were not pre-allocated
_pushgrid!(::Type{<:NamedTuple}, o) = push!(o, map(grid -> similar(grid), o[1]))
_pushgrid!(::Type{<:AbstractArray}, o) = push!(o, similar(o[1]))

function initialise!(o::GraphicOutput, data) 
    initalisegraphics(o, data)
end
function finalise!(o::GraphicOutput, data) 
    _storeframe!(o, eltype(o), data)
    finalisegraphics(o, data)
end

# Additional interface for GraphicOutput

function showframe(o::GraphicOutput, data) 
    showframe(o, proc(data), data)
end
function showframe(o::GraphicOutput, ::Processor, data)
    frame = map(maybe_rebuild_grid, init(o), grids(data))
    showframe(frame, o, data)
end
# Handle NamedTuple for outputs that only accept AbstractArray
showframe(frame::NamedTuple, o::GraphicOutput, data) = showframe(first(frame), o, data)

function maybe_rebuild_grid(a, b)
    if a isa AbstractDimArray 
        b isa AbstractDimArray ? rebuild(a, parent(b)) : rebuild(a, b) 
    else
        b
    end
end

initalisegraphics(o::GraphicOutput, data) = nothing
finalisegraphics(o::GraphicOutput, data) = showframe(o, data)
