"""
    Output

Abstract supertype for simulation outputs.

Outputs are store or display simulation results, usually
as a vector of grids, one for each timestep - but they may also
sum, combine or otherwise manipulate the simulation grids to improve
performance, reduce memory overheads or similar.

Simulation outputs are decoupled from simulation behaviour,
and in many cases can be used interchangeably.
"""
abstract type Output{T,A} <: AbstractDimArray{T,1,Tuple{Ti},A} end
# Generic ImageOutput constructor. Converts an init array to vector of arrays.
function (::Type{T})(
    init::Union{NamedTuple,AbstractArray}; extent=nothing, kw...
) where T <: Output
    extent = extent isa Nothing ? Extent(; init=init, kw...) : extent
    frames = [_replicate_init(init, replicates(extent))]
    T(; frames, running=false, extent=extent, kw...)
end

function Base.summary(io::IO, A::Output)
    DimensionalData.print_ndims(io, size(A))
    print(io, string(nameof(typeof(A))))
    DimensionalData.print_name(io, name(A))
end

function DimensionalData.show_after(io::IO, mime, A::Output)
    blockwidth = get(io, :blockwidth, 0)
    DimensionalData.print_block_close(io, blockwidth)
    ndims(A) > 0 && println(io)
    if eltype(A) <: NamedTuple
        simplified_grids = map(A) do f
            DimensionalData.ShowWith(join(map((k, v) -> "$k=$(nameof(typeof(v)))", keys(f), values(f)), ", "))
        end
        DimensionalData.print_array(io, mime, simplified_grids)
    else
        DimensionalData.print_array(io, mime, A)
    end
end

# Forward base methods to the frames array
Base.parent(o::Output) = frames(o)
Base.length(o::Output) = length(parent(o))
Base.size(o::Output) = size(parent(o))
Base.firstindex(o::Output) = firstindex(parent(o))
Base.lastindex(o::Output) = lastindex(parent(o))
Base.push!(o::Output, x) = push!(parent(o), x)
Base.step(o::Output) = step(tspan(o))

# DimensionalData interface ######################################################
# This allows indexing the output using values from tspan
function DimensionalData.dims(o::Output)
    ts = tspan(o)
    val = isstored(o) ? ts : ts[end]:step(ts):ts[end]
    (Ti(Sampled(val; order=ForwardOrdered(), span=Regular(step(ts)), sampling=Intervals(Start()))),)
end
DimensionalData.refdims(o::Output) = ()
DimensionalData.name(o::Output) = Symbol("")
DimensionalData.metadata(o::Output) = NoMetadata()
# Output rebuild just returns a DimArray
function DimensionalData.rebuild(o::Output, data, dims::Tuple, refdims, name, metadata)
    return DimArray(data, dims, refdims, name, metadata) 
end
function DimensionalData.rebuild(o::Output; 
    data=frames(o), dims=DD.dims(o), refdims=DD.refdims(o), names=DD.name(o), metadata=DD.metadata(o), kw...
)
    return DimArray(data, dims, refdims, name, metadata) 
end

# Required getters and setters for all outputs ###################################
frames(o::Output) = o.frames
isrunning(o::Output) = o.running
extent(o::Output) = o.extent
init(o::Output) = init(extent(o))
mask(o::Output) = mask(extent(o))
aux(o::Output, key...) = aux(extent(o), key...)
padval(o::Output, key...) = padval(extent(o), key...)
replicates(o::Output) = replicates(extent(o))
tspan(o::Output) = tspan(extent(o))
timestep(o::Output) = step(tspan(o))

setrunning!(o::Output, val) = o.running = val
settspan!(o::Output, tspan) = settspan!(extent(o), tspan)

# Default values for optional getters and setters ################################
ruleset(o::Output) =
    throw(ArgumentError("No ruleset on the output. Pass one to `sim!` as the second argument"))
fps(o::Output) = nothing
stoppedframe(o::Output) = lastindex(o)

setfps!(o::Output, x) = nothing
settimestamp!(o::Output, f) = nothing
setstoppedframe!(o::Output, f) = nothing

isasync(o::Output) = false
isstored(o::Output) = true
isshowable(o::Output, frame) = false
initialise!(o::Output, data) = nothing
finalise!(o::Output, data) = nothing
initialisegraphics(o::Output, data) = nothing
finalisegraphics(o::Output, data) = nothing
maybesleep(o::Output, frame) = sleep(0) # Helps for exiting simulations
showframe(o::Output, data) = nothing

frameindex(o::Output, data::AbstractSimData) = frameindex(o, currentframe(data))
frameindex(o::Output, f::Int) = isstored(o) ? f : oneunit(f)


# Storing grid values to outputs during the simulation ###################################

function storeframe!(output::Output, data)
    # Make sure the frame exists in the output
    checkbounds(output, frameindex(output, data))
    # copy to it
    _storeframe!(output, eltype(output), data)
end

# copy one or multiple frames from grid/s to the Output
function _storeframe!(output::Output, ::Type{<:AbstractArray}, data)
    grid = first(grids(data))
    _copyto_output!(output[frameindex(output, data)], grid)
end
function _storeframe!(output::Output, ::Type{<:NamedTuple}, data)
    map(values(grids(data)), keys(data)) do grid, key
        _copyto_output!(output[frameindex(output, data)][key], grid)
    end
end

# Copy cells from grid to output
_copyto_output!(outgrid, grid::GridData) = _copyto_output!(outgrid, grid, proc(grid))
function _copyto_output!(outgrid, grid::GridData, proc::CPU)
    copyto!(outgrid, CartesianIndices(outgrid), grid, CartesianIndices(outgrid))
end
# Copy cells from grid to output using multiple threads
function _copyto_output!(outgrid, grid::GridData{<:Any,Tuple{X,Y}}, proc::ThreadedCPU) where {X, Y}
    Threads.@threads for j in axes(outgrid, ndims(outgrid)) 
        for i in axes(outgrid, 1)
            @inbounds outgrid[i, j] = grid[i, j]
        end
    end
end

# Grid initialisation ###################################################################

# Grids are preallocated and reused.
init_output_grids!(o::Output, init) = init_output_grids!(o[1], o::Output, init)
# Array grids are copied
function init_output_grids!(grid::AbstractArray, o::Output, init::AbstractArray)
    copyto!(grid, init)
    return o
end
# All arrays are copied if both are named tuples
function init_output_grids!(grids::NamedTuple, o::Output, inits::NamedTuple)
    map(grids, inits) do grid, init
        copyto!(grid, init)
    end
    return o
end
