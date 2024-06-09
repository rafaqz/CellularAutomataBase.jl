
"""
    SparseOpt <: PerformanceOpt

    SparseOpt()

An optimisation flag that ignores all padding valuesin the grid,
by default zeros.

For low-density simulations performance may improve by
orders of magnitude, as only used cells are run.

Specifiy with:

```julia
ruleset = Ruleset(rule; opt=SparseOpt())
# or
output = sim!(output, rule; opt=SparseOpt())
```

`SparseOpt` is best demonstrated with this simulation, where the grey areas do not
run except where the stencil partially hangs over an area that is not grey:

![SparseOpt demonstration](https://raw.githubusercontent.com/cesaraustralia/DynamicGrids.jl/media/complexlife_spareseopt.gif)
"""
struct SparseOpt{F<:Function} <: PerformanceOpt
    f::F
end
SparseOpt() = SparseOpt(==(0))

@inline _isactive(val, opt::SparseOpt{<:Function}) = !opt.f(val)

sourcestatus(d::GridData) = optdata(d).sourcestatus
deststatus(d::GridData) = optdata(d).deststatus

# Run kernels with SparseOpt, block by block:
function map_with_optimisation(
    f, simdata::AbstractSimData{S}, proc, ::SparseOpt, ruletype::Val{<:Rule}, rkeys
) where {S<:Tuple{Y,X}} where {Y,X}
    # Only use SparseOpt for single-grid rules with grid radii > 0
    grid = _firstgrid(simdata, rkeys)
    R = radius(grid)
    if R == 0
        map_with_optimisation(f, simdata, proc, NoOpt(), ruletype, rkeys)
        return nothing
    end
    B = 2R
    status = sourcestatus(grid)
    let f=f, proc=proc, status=status
        map_on_processor(proc, simdata, 1:_indtoblock(X+R, B)) do bj
            for  bi in 1:_indtoblock(Y+R, B)
                status[bi, bj] || continue
                # Convert from padded block to init dimensions
                istart, jstart = _blocktoind(bi, B) - R, _blocktoind(bj, B) - R
                # Stop at the init row/column size, not the padding or block multiple
                istop, jstop = min(istart + B - 1, Y), min(jstart + B - 1, X)
                # Skip the padding
                istart, jstart  = max(istart, 1), max(jstart, 1)
                for j in jstart:jstop 
                    @simd for i in istart:istop
                        f((i, j))
                    end
                end
            end
        end
    end
    return nothing
end

function row_kernel!(
    simdata::AbstractSimData, grid::GridData{<:Any,<:Tuple{Y,X},R}, proc, opt::SparseOpt,
    ruletype::Val, rule::Rule, rkeys, wkeys, bi
) where {Y,X,R}
    # No SparseOpt for radius 0
    if R === 0
        return row_kernel!(simdata, grid, proc, NoOpt(), ruletype, rule, rkeys, wkeys, bi)
    end
    B = 2R
    S = 2R + 1
    nblockcols = _indtoblock(X+R, B)
    srcstatus, dststatus = sourcestatus(grid), deststatus(grid)

    # Blocks ignore padding! the first block contains padding.
    i = _blocktoind(bi, B)
    i > Y && return nothing
    # Get current bloc
    skippedlastblock = true

    # Initialise block status for the start of the row
    # The first column always runs, it's buggy otherwise.
    @inbounds bs11, bs12 = true, true
    @inbounds bs21, bs22 = true, true
    # New block status
    newbs12 = false
    newbs22 = false
    for bj = 1:nblockcols
        # Shuffle current window status
        bs11, bs21 = bs12, bs22
        @inbounds bs12, bs22 = srcstatus[bi, bj + 1], srcstatus[bi + 1, bj + 1]
        # Skip this block if it and the neighboring blocks are inactive
        if !(bs11 | bs12 | bs21 | bs22)
            skippedlastblock = true
            # Run the rest of the chain if it exists and more than 1 grid is used
            if rule isa Chain && length(rule) > 1 && length(rkeys) > 1
                # Loop over the grid COLUMNS inside the block
                jstart = _blocktoind(bj, B)
                jstop = min(jstart + B - 1, X)
                for j in jstart:jstop
                    # Loop over the grid ROWS inside the block
                    blocklen = min(Y, i + B - 1) - i + 1
                    for b in 1:blocklen
                        cell_kernel!(simdata, ruletype, rule, rkeys, wkeys, i + b - 1, j)
                    end
                end
            end
            continue
        end
        # Define area to loop over with the block.
        # It's variable because the last block may be partial
        jstart = _blocktoind(bj, B)
        jstop = min(jstart + B - 1, X)

        # Reinitialise stencil windows if we have skipped a section of the array
        if skippedlastblock
            skippedlastblock = false
        end
        # Shuffle new window status
        newbs11 = newbs12
        newbs21 = newbs22
        newbs12 = newbs22 = false

        # Loop over the grid COLUMNS inside the block
        for j in jstart:jstop
            # Update windows unless feshly populated
            # Which block column are we in, 1 or 2
            curblockj = (j - jstart) ÷ R + 1
            # Loop over the COLUMN of windows covering the block
            blocklen = min(Y, i + B - 1) - i + 1
            for b in 1:blocklen
                # Set rule window
                rule1 = Stencils.rebuild(rule, unsafe_neighbors(stencil(rule), grid, CartesianIndex(i, j)))
                # Run the rule kernel for the cell
                writeval = cell_kernel!(simdata, ruletype, rule1, rkeys, wkeys, i + b - 1, j)
                # Update the status for the current block
                cs = _cellstatus(opt, wkeys, writeval)
                curblocki = R == 1 ? b : (b - 1) ÷ R + 1
                if curblocki == 1
                    curblockj == 1 ? (newbs11 |= cs) : (newbs12 |= cs)
                else
                    curblockj == 1 ? (newbs21 |= cs) : (newbs22 |= cs)
                end
            end

            # Combine new block status with deststatus array
            @inbounds dststatus[bi, bj] |= newbs11
            @inbounds dststatus[bi+1, bj] |= newbs21
            @inbounds dststatus[bi, bj+1] |= newbs12
            @inbounds dststatus[bi+1, bj+1] |= newbs22
        end
    end
    return nothing
end

@inline _cellstatus(opt::SparseOpt, wkeys::Tuple, writeval) = _isactive(writeval[1], opt)
@inline _cellstatus(opt::SparseOpt, wkeys, writeval) = _isactive(writeval, opt)

# SparseOpt methods

# _build_status => (Matrix{Bool}, Matrix{Bool})
# Build block-status arrays
# We add an additional block that is never used so we can 
# index into it in the block loop without checking
function _build_optdata(opt::SparseOpt, source, r::Int)
    r > 0 || return nothing
    hoodsize = 2r + 1
    blocksize = 2r
    nblocs = _indtoblock.(size(source), blocksize) .+ 1
    sourcestatus = zeros(Bool, nblocs)
    deststatus = zeros(Bool, nblocs)
    return (; sourcestatus, deststatus)
end

Stencils.switch(::SparseOpt, ::Nothing) = nothing
function Stencils.switch(::SparseOpt, optdata)
    (sourcestatus=optdata.deststatus, deststatus=optdata.sourcestatus)
end


# Initialise the block status array.
# This tracks whether anything has to be done in an area of the main array.
function _update_optdata!(grid, opt::SparseOpt)
    isnothing(optdata(grid)) && return grid
    blocksize = 2 * radius(grid)
    src = parent(source(grid))
    for I in CartesianIndices(src)
        # Mark the status block (by default a non-zero value)
        if _isactive(src[I], opt)
            bi = _indtoblock.(Tuple(I), blocksize)
            @inbounds sourcestatus(grid)[bi...] = true
            @inbounds deststatus(grid)[bi...] = true
        end
    end
    return grid
end

# Clear the destination grid and its status for SparseOpt.
# NoOpt writes to every cell, so this is not required
function _cleardest!(grid, opt::SparseOpt)
    dest(grid) .= source(grid)
    if !isnothing(optdata(grid))
        deststatus(grid) .= false
    end
end

# Sets the status of the destination block that the current index is in.
# It can't turn of block status as the block is larger than the cell
# But should be used inside a LOCK
function _setoptindex!(grid::GridData{<:SwitchMode,<:Any,R}, opt::SparseOpt, x, I...) where R
    isnothing(optdata(grid)) && return grid
    blockindex = _indtoblock.(I .+ R, 2R)
    @inbounds deststatus(grid)[blockindex...] |= !(opt.f(x))
    return nothing
end

# _wrapopt!
# Copies status from opposite sides/corners in Wrap boundary mode
function Stencils.after_update_boundary!(grid::GridData, ::SparseOpt)
    isnothing(optdata(grid)) && return grid

    status = sourcestatus(grid)
    # !!! The end row/column is always empty !!!
    # Its padding for block opt. So we work with end-1 and end-2
    # We should probably re-write this using the known grid sizes,
    # instead of `end`
    
    # This could be further optimised by not copying the end-2 
    # block column/row when the blocks are aligned at both ends.
    
    # Sides
    status[1, :] .|= status[end-1, :] .| status[end-2, :]
    status[:, 1] .|= status[:, end-1] .| status[:, end-2]
    status[end-1, :] .|= status[1, :]
    status[:, end-1] .|= status[:, 1]
    status[end-2, :] .|= status[1, :]
    status[:, end-2] .|= status[:, 1]

    # Corners
    status[1, 1] |= status[end-1, end-1] | status[end-2, end-1] | 
                    status[end-1, end-2] | status[end-2, end-2]
    status[end-1, 1] |= status[1, end-1] | status[1, end-2]
    status[end-2, 1] |= status[1, end-1] | status[1, end-2]
    status[1, end-1] |= status[end-1, 1] | status[end-2, 1]
    status[1, end-2] |= status[end-1, 1] | status[end-2, 1]
    status[end-1, end-1] |= status[1, 1] 
    status[end-2, end-2] |= status[1, 1] 

    return grid
end


"""
    SparseOptInspector()

A [`Renderer`](@ref) that checks [`SparseOpt`](@ref) visually.
Cells that do not run show in gray. Errors show in red, but if they do there's a bug.
"""
struct SparseOptInspector{A} <: SingleGridRenderer
    accessor::A
end
SparseOptInspector() = SparseOptInspector(identity)

accessor(p::SparseOptInspector) = p.accessor

function cell_to_pixel(p::SparseOptInspector, mask, minval, maxval, data::AbstractSimData, val, I::Tuple)
    opt(data) isa SparseOpt || error("Can only use SparseOptInspector with SparseOpt grids")
    r = radius(first(grids(data)))
    blocksize = 2r
    blockindex = _indtoblock.((I[1] + r,  I[2] + r), blocksize)
    normedval = normalise(val, minval, maxval)
    # This is done at the start of the next frame, so wont show up in
    # the image properly. So do it preemtively?
    Stencils.after_update_boundary!(first(data))
    status = sourcestatus(first(data))
    if status[blockindex...]
        if normedval > 0
            to_rgb(normedval)
        else
            to_rgb((0.0, 0.0, 0.0))
        end
    elseif normedval > 0
        to_rgb((1.0, 0.0, 0.0)) # This (a red cell) would mean there is a bug in SparseOpt
    else
        to_rgb((0.5, 0.5, 0.5))
    end
end

# Custom SimSettings constructor: SparseOpt does not work on GPU
function SimSettings(
    boundary::B, proc::P, opt::SparseOpt, cellsize::C, timestep::T
) where {B,P<:GPU,C,T}
    @info "SparseOpt does not work on GPU. Using NoOpt instead."
    SimSettings{B,P,NoOpt,C,T}(boundary, proc, NoOpt(), cellsize, timestep)
end
