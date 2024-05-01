abstract type GridMode end
abstract type ReadMode <: GridMode end
abstract type WriteMode <: GridMode end
abstract type SwitchMode <: WriteMode end

"""
    AbstractGridArray <: AbstractSwitchingStencilArray

Simulation data specific to a single grid.

These behave like arrays, but contain both source and
destination arrays as simulations need separate read and
write steps to maintain independence between cells.

`GridData` objects also contain other data and settings needed
for optimisations.

# Type parameters

- `S`: grid size type tuple
- `R`: grid padding radius
- `T`: grid data type
"""
abstract type AbstractGridData{Mode,S,R,T,N,A,H,B,P} <: Stencils.AbstractSwitchingStencilArray{S,R,T,N,A,H,B,P} end

# Getters
proc(d::AbstractGridData) = d.proc
opt(d::AbstractGridData) = d.opt
optdata(d::AbstractGridData) = d.optdata
mask(d::AbstractGridData) = d.mask
maskval(d::AbstractGridData) = d.maskval
replicates(d::AbstractGridData) = d.replicates
indices(d::AbstractGridData) = d.indices

# Get the size of the grid
gridsize(d::AbstractGridData) = size(d)
gridsize(A::AbstractArray) = size(A)
gridsize(nt::NamedTuple) = gridsize(first(nt))
gridsize(nt::NamedTuple{(),Tuple{}}) = 0, 0

# Get a view of the grid, without padding
gridview(d::AbstractGridData) = sourceview(d)
# Get a view of the grid source, without padding
sourceview(d::AbstractGridData) = _unpad_view(source(d), d)
# Get a view of the grid dest, without padding
destview(d::AbstractGridData) = _unpad_view(dest(d), d)

_unpad_view(A, d::AbstractGridData{<:Any,<:Any,R}) where R =
    view(A, map(a -> a .+ R, axes(d))...)

# Get an a view of the source, preferring the underlying array if there is no radius
source_array_or_view(d::AbstractGridData) = sourceview(d)
source_array_or_view(d::AbstractGridData{<:Any,<:Any,0}) = source(d)
# Get an a view of the dest, preferring the underlying array if there is no radius
dest_array_or_view(d::AbstractGridData) = destview(d)
dest_array_or_view(d::AbstractGridData{<:Any,<:Any,0}) = dest(d)

_build_optdata(opt::PerformanceOpt, init, r) = nothing

# _indtoblock
# Convert regular index to block index
@inline _indtoblock(x::Int, blocksize::Int) = (x - 1) ÷ blocksize + 1

# _blocktoind
# Convert block index to regular index
@inline _blocktoind(x::Int, blocksize::Int) = (x - 1) * blocksize + 1

Base.copy!(dest::AbstractGridData, source::AbstractDimArray) = copy!(dest, parent(source))

"""
    GridData <: GridData

    GridData(grid::GridData)
    GridData{S,R}(init::AbstractArray, mask, opt, boundary, padding)

[`GridData`](@ref) object passed to rules for reading only.

Has `ReadMode`, `WriteMode` and `SwitchMode` to control behaviour.

Reads are always from the `source` array.
"""
struct GridData{
    Mode,S<:Tuple,R,T,N,A,H,B,P,Pr<:Processor,Op<:PerformanceOpt,OpD,Ma,MV,Re,I
} <: AbstractGridData{Mode,S,R,T,N,A,H,B,P}
    source::A
    dest::A
    stencil::H
    boundary::B
    padding::P
    proc::Pr
    opt::Op
    optdata::OpD
    mask::Ma
    maskval::MV
    replicates::Re
    indices::I
end
function GridData{Mode,S,R}(
    source::A, dest::A, stencil::H, boundary::B, padding::P,
    proc::Pr, opt::Op, optdata::OpD, mask::Ma, maskval::MV,
    replicates::Re, indices::I
) where {Mode,S,R,A<:AbstractArray{T,N},H,B,P,Pr,Op,OpD,Ma,MV,Re,I} where {T,N}
    GridData{Mode,S,R,T,N,A,H,B,P,Pr,Op,OpD,Ma,MV,Re,I}(
        source, dest, stencil, boundary, padding, proc, opt, optdata, mask, maskval, replicates, indices
    )
end
function GridData{Mode}(d::AbstractGridData{<:Any,S,R,T,N,A}) where {Mode,S,R,T,N,A}
    args = stencil(d), boundary(d), padding(d), proc(d), opt(d), optdata(d), mask(d), maskval(d), replicates(d), indices(d)
    GridData{Mode,S,R,T,N,A,map(typeof, args)...}(source(d), dest(d), args...)
end
@inline function GridData{Mode,S,R}(
    init::AbstractArray{<:Any,N}, stencil::Stencil, boundary::BoundaryCondition,
    padding::Padding, proc, opt, mask, maskval, replicates, indices=nothing
) where {Mode,S,R,N}
    # If the grid radius is larger than zero we pad it
    if R > 0
        # Blocks (only used for 2d sims) need additional vertical padding.
        # TODO: this needs clarification.
        r = N == 2 ? ((R, 3R), (R, R)) : R
        source = Stencils.pad_array(padding, boundary, stencil, init)
        dest = Stencils.pad_array(padding, boundary, stencil, init)
        # TODO is there a chance this doen't make a copy?
        # For now just assert that they are not identical arrays
        @assert source !== dest
    else
        source = deepcopy(init)
        dest = deepcopy(init)
    end
    optdata = _build_optdata(opt, source, R)

    grid = GridData{Mode,S,R}(
        source, dest, stencil, boundary, padding, proc, opt, optdata, mask, maskval, replicates, indices
    )
    update_boundary!(grid)
    return grid
end

function ConstructionBase.constructorof(
    ::Type{G}
) where G<:GridData{Mode,S,R} where {Mode,S,R}
    GridData{Mode,S,R}
end

for f in (:getindex, :view, :dotview)
    unsafe_f = Symbol(string("unsafe_", f))
    @eval begin
        Base.@propagate_inbounds function Base.$f(A::GridData, i1::Int, Is::Int...)
            # If we have replicates, we add the replicate index here
            I = _maybe_complete_indices(A, (i1, Is...))
            Stencils.$unsafe_f(A, I...)
        end
    end
end

### UNSAFE / LOCKS required

# Base.setindex!
# This is not safe for general use.
# It can be used where only identical transformations of a cell
# can happen from any other cell, such as setting all 1s to 2.
@propagate_inbounds function Base.setindex!(d::GridData{<:WriteMode}, x, i1::Int, Is::Int...)
    I = _maybe_complete_indices(d, add_halo(d, (i1, Is...)))
    _setindex!(d, proc(d), x, I...)
end

@propagate_inbounds function _setindex!(d::GridData{<:WriteMode}, proc::SingleCPU, x, I...)
    @boundscheck checkbounds(dest(d), I...)
    # @inbounds
    _setoptindex!(d, x, I...)
    source(d)[d...] = x
end
@propagate_inbounds function _setindex!(d::GridData{<:WriteMode}, proc::ThreadedCPU, x, I...)
    # Dest status is not threadsafe, even if the
    # setindex itself is safe. So we LOCK
    lock(proc)
    # @inbounds
    _setoptindex!(d, x, I...)
    unlock(proc)
    source(d)[d...] = x
end
@propagate_inbounds function _setindex!(d::GridData{<:SwitchMode}, proc::SingleCPU, x, I...)
    @boundscheck checkbounds(dest(d), I...)
    # @inbounds
    _setoptindex!(d, x, I...)
    dest(d)[I...] = x
end
@propagate_inbounds function _setindex!(d::GridData{<:SwitchMode}, proc::ThreadedCPU, x, I...)
    # Dest status is not threadsafe, even if the
    # setindex itself is safe. So we LOCK
    lock(proc)
    # @inbounds
    _setoptindex!(d, x, I...)
    unlock(proc)
    dest(d)[I...] = x
end

function _maybe_complete_indices(data::GridData, I::Tuple)
    if isnothing(replicates(data)) || isnothing(indices(data))
        I
    else
        (I..., last(indices(data)))
    end
end

# _setoptindex!
# Do anything the optimisation needs on `setindex`
_setoptindex!(d::GridData{<:WriteMode,<:Any,R}, x, I...) where R = _setoptindex!(d, opt(d), x, I...)
_setoptindex!(d::GridData{<:WriteMode,<:Any,R}, opt::PerformanceOpt, x, I...) where R = nothing

function Stencils.switch(grids::NamedTuple{<:Any,Tuple{T,Vararg}}) where {T<:GridData{<:SwitchMode}}
    map(switch, grids)
end
function Stencils.switch(A::T) where {T<:GridData{<:SwitchMode}}
    od = switch(opt(A), optdata(A))
    T(dest(A), source(A), stencil(A), boundary(A), padding(A), proc(A), opt(A), od, mask(A), maskval(A), replicates(A), indices(A))
end
Stencils.switch(::PerformanceOpt, optdata) = optdata

Stencils.after_update_boundary!(grid::GridData) = Stencils.after_update_boundary!(grid, opt(grid))
Stencils.after_update_boundary!(grid::GridData, opt) = grid

function Base.copy!(S::GridData{<:Any,R}, A::AbstractDimArray) where R
    pad_axes = add_halo(S, axes(S))
    copyto!(parent(parent(S)), CartesianIndices(pad_axes), A, CartesianIndices(A))
    return
end
function Base.copy!(A::AbstractDimArray, S::GridData{<:Any,R}) where R
    pad_axes = add_halo(S, axes(S))
    copyto!(A, CartesianIndices(A), parent(S), CartesianIndices(pad_axes))
    return A
end
function Base.copy!(dst::GridData{<:Any,RD}, src::GridData{<:Any,RS}) where {RD,RS}
    dst_axes = add_halo(dst, axes(dst))
    src_axes = add_halo(src, axes(src))
    copyto!(
        parent(dst), CartesianIndices(dst_axes),
        parent(src), CartesianIndices(src_axes)
    )
    return dst
end

function Base.copyto!(dst::GridData, idst::CartesianIndices, src::GridData, isrc::CartesianIndices)
    dst_axes = add_halo(dst, idst)
    src_axes = add_halo(src, isrc)
    copyto!(
        parent(dst), CartesianIndices(dst_axes),
        parent(src), CartesianIndices(src_axes),
    )
    return dst
end
function Base.copyto!(dst::AbstractDimArray, idst::CartesianIndices, src::GridData, isrc::CartesianIndices)
    src_axes = add_halo(src, isrc)
    copyto!(dst, idst, parent(src), src_axes)
    return dst
end
function Base.copyto!(dst::GridData, idst::CartesianIndices, src::AbstractDimArray, isrc::CartesianIndices)
    dst_axes = add_halo(dst, idst)
    copyto!(source(dst), dst_axes, src, isrc)
    return dst
end
function Base.copyto!(dst::GridData, src::GridData)
    dst_axes = add_halo(dst, axes(dst))
    src_axes = add_halo(src, axes(src))
    copyto!(
        source(dst), CartesianIndices(dst_axes),
        source(src), CartesianIndices(src_axes),
    )
    return dst
end
function Base.copyto!(dst::AbstractDimArray, src::GridData)
    src_axes = add_halo(src, axes(src))
    copyto!(
        dst, CartesianIndices(dst),
        source(src), CartesianIndices(src_axes),
    )
    return dst
end
function Base.copyto!(dst::GridData, src::AbstractDimArray)
    dst_axes = add_halo(dst, axes(dst))
    copyto!(
        source(dst), CartesianIndices(dst_axes),
        src, CartesianIndices(src),
    )
    return dst
end
