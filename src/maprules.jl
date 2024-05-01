
# Map a rule over the grids it reads from, updating the grids it writes to.
#
# This is split into setup methods and application methods,
# for dispatch and to introduce a function barrier for type stability.

# We dispatch on `ruletype(rule)` to allow wrapper rules
# to pass through the type of the wrapped rule.
# Putting the type in `Val` is best for performance.
maprule!(data::AbstractSimData, rule) =
    maprule!(data, _val_ruletype(rule), rule)

# CellRule
function maprule!(simdata::AbstractSimData, ruletype::Val{<:CellRule}, rule)
    rkeys, _ = _getreadgrids(rule, simdata)
    wkeys, _ = _getwritegrids(WriteMode, rule, simdata)
    ruledata = RuleData(simdata, rule)
    maprule!(ruledata, ruletype, rule, rkeys, wkeys)
    # Not swictch required
    return simdata
end

# NeighborhoodRule
function maprule!(simdata::AbstractSimData, ruletype::Val{<:NeighborhoodRule}, rule)
    rkeys, rgrids = _getreadgrids(rule, simdata)
    wkeys, wgrids = _getwritegrids(SwitchMode, rule, simdata)
    # Copy or zero out boundary where needed
    _update_boundary!(rgrids)
    _cleardest!(simdata[stencilkey(rule)])
    ruledata = RuleData(simdata, rule)
    maprule!(ruledata, ruletype, rule, rkeys, wkeys)
    # Swap the dest/source of grids that were written to
    # and combine the written grids with the original simdata
    final_rgrids = _to_readonly(_switch(wgrids))
    return _replacegrids(simdata, wkeys, final_rgrids)
end

# SetRule
function maprule!(simdata::AbstractSimData, ruletype::Val{<:SetRule}, rule)
    rkeys, _ = _getreadgrids(rule, simdata)
    wkeys, wgrids = _getwritegrids(SwitchMode, rule, simdata)
    map(_astuple(wkeys, wgrids)) do g
        copyto!(parent(dest(g)), parent(source(g)))
    end
    ruledata = RuleData(_combinegrids(simdata, wkeys, wgrids), rule)
    maprule!(ruledata, ruletype, rule, rkeys, wkeys)
    # We have to mask after SetRule
    final_grids = _maybemask!(_to_readonly(_switch(wgrids)))
    return _replacegrids(simdata, wkeys, final_grids)
end

_switch(xs::Tuple) = map(switch, xs)
_switch(x) = switch(x)

# SetGridRule
function maprule!(simdata::AbstractSimData, ruletype::Val{<:SetGridRule}, rule)
    rkeys, rgrids = _getreadgrids(rule, simdata)
    wkeys, wgrids = _getwritegrids(WriteMode, rule, simdata)
    ruledata = RuleData(_combinegrids(simdata, wkeys, wgrids), rule)
    # Run the rule
    applyrule!(ruledata, rule)
    # We don't mask here or do anything, its on the user
    return simdata
end

# Expand method arguments for dispatch on processor and optimisation
function maprule!(ruledata::RuleData, ruletype::Val, rule, rkeys, wkeys)
    maprule!(ruledata, proc(ruledata), opt(ruledata), ruletype, rule, rkeys, wkeys)
end

_update_boundary!(gs::Union{NamedTuple,Tuple}) = map(_update_boundary!, gs)
_update_boundary!(g::GridData) = update_boundary!(g)

# Most Rules
# 2 dimensional, with processor selection and optimisations in `map_with_optimisation`
function maprule!(ruledata::RuleData{<:Any,<:Tuple{Y,X}}, proc::CPU, opt, ruletype::Val, rule, rkeys, wkeys) where {Y,X}
    let ruledata=ruledata, proc=proc, opt=opt, rule=rule,
        rkeys=rkeys, wkeys=wkeys, ruletype=ruletype
        map_with_optimisation(ruledata, proc, opt, ruletype, rkeys) do I 
            cell_kernel!(data, ruletype, rule, rkeys, wkeys, I...)
        end
    end
end
# Arbitrary dimensions, no processor or optimisation selection beyond CPU/GPU
function maprule!(ruledata::RuleData, proc::CPU, opt, ruletype::Val, rule, rkeys, wkeys)
    let ruledata=ruledata, proc=proc, opt=opt, rule=rule,
        rkeys=rkeys, wkeys=wkeys, ruletype=ruletype
        for I in CartesianIndices(first(grids(ruledata)))
            cell_kernel!(ruledata, ruletype, rule, rkeys, wkeys, Tuple(I)...)
        end
    end
end

# Neighborhood rules
# 2 dimensional, with processor selection and optimisations
# function maprule!(
#     data::AbstractSimData{<:Any,<:Tuple{I,J}}, proc::CPU, opt, ruletype::Val{<:NeighborhoodRule}, 
#     rule, rkeys, wkeys
# ) where {I,J}
#     hoodgrid = _firstgrid(data, rkeys)
#     let data=data, hoodgrid=hoodgrid, proc=proc, opt=opt, ruletyp=ruletype, rule=rule, rkeys=rkeys, wkeys=wkeys
#         B = 2radius(hoodgrid)
#         # UNSAFE: we must avoid sharing status blocks, it could cause race conditions 
#         # when setting status from different threads. So we split the grid in 2 interleaved
#         # sets of rows, so that we never run adjacent rows simultaneously
#         map_on_processor(proc, data, 1:2:_indtoblock(I, B)) do bi
#             row_kernel!(data, hoodgrid, proc, opt, ruletype, rule, rkeys, wkeys, bi)
#         end
#         map_on_processor(proc, data, 2:2:_indtoblock(I, B)) do bi
#             row_kernel!(data, hoodgrid, proc, opt, ruletype, rule, rkeys, wkeys, bi)
#         end
#     end
#     return nothing
# end
# Arbitrary dimensions, no processor or optimisation selection beyond CPU/GPU
function maprule!(
    data::RuleData, proc::SingleCPU, opt, ruletype::Val{<:NeighborhoodRule}, rule, rkeys, wkeys
)
    hoodgrid = _firstgrid(data, rkeys)
    for I in CartesianIndices(hoodgrid) 
        stencil_kernel!(data, hoodgrid, ruletype, rule, rkeys, wkeys, Tuple(I)...)
    end
    return nothing
end
function maprule!(
    data::RuleData, proc::ThreadedCPU, opt, ruletype::Val{<:NeighborhoodRule}, rule, rkeys, wkeys
)
    hoodgrid = _firstgrid(data, rkeys)
    map_with_optimisation(data, proc, opt, ruletype, rkeys) do I
        stencil_kernel!(data, hoodgrid, ruletype, rule, rkeys, wkeys, Tuple(I)...)
    end
    return nothing
end

### Rules that don't need a stencil window ####################

# map_with_optimisation
# Map kernel over the grid, specialising on PerformanceOpt.

# Run kernel over the whole grid, cell by cell:
function map_with_optimisation(
    f, simdata::AbstractSimData{S}, proc, ::NoOpt, ::Val{<:Rule}, rkeys
) where S<:Tuple{I,J} where {I,J}
    map_on_processor(proc, simdata, 1:J) do j
        @simd for i in 1:I
            f((i, j))
        end
    end
end
function map_with_optimisation(
    f, simdata::AbstractSimData{S}, proc, ::NoOpt, ::Val{<:Rule}, rkeys
) where S<:Tuple{I,J,K} where {I,J,K}
    map_on_processor(proc, simdata, 1:K) do k
        for j in 1:J 
            @simd for i in 1:I
                f((i, j, k))
            end
        end
    end
end
function map_with_optimisation(
    f, simdata::AbstractSimData{S}, proc, ::NoOpt, ::Val{<:Rule}, rkeys
) where S<:Tuple{I,J,K,L} where {I,J,K,L}
    map_on_processor(proc, simdata, 1:L) do l
        for k in 1:K 
            for j in 1:J 
                @simd for i in 1:I
                    f((i, j, k, l))
                end
            end
        end
    end
end

# broacast_on_processor
# Map kernel over the grid, specialising on the processor
#
# Looping over cells or blocks on a single CPU
@inline function map_on_processor(f, proc::SingleCPU, data, range)
    for n in range
        f(n) # Run rule over each column
    end
end
# Or threaded on multiple CPUs
@inline function map_on_processor(f, proc::ThreadedCPU, data, range)
    Threads.@threads :static for n in range
        f(n) # Run rule over each column
    end
    # We don't want to share memory between neighborhoods
    # min_cols = max(3, 2radius(data) + 1)
    # N = Threads.nthreads()
    # allchunks = collect(Iterators.partition(rnge, min_cols))
    # chunks = map(1:N) do i
    #     allchunks[i:N:end]
    # end
    # tasks = map(chunks) do chunk
    #     Threads.@spawn begin
    #         for subchunk in chunk
    #             for n in subchunk
    #                 f(n)
    #             end
    #         end
    #     end
    # end
    # states = fetch.(tasks)
    return nothing
end

# cell_kernel!
# runs a rule for the current cell
@inline function cell_kernel!(data::RuleData, ruletype, rule, rkeys, wkeys, I...)
    # When we have replicates as an additional grid 
    # dimension we hide the extra dimension from rules.
    I1 = _strip_replicates(data, I)
    # We skip the cell if there is a mask layer
    m = mask(data)
    if !isnothing(m)
        m[I1...] || return nothing
    end
    # We read a value from the grid
    readval = _readcell(data, rkeys, I...)
    # Update the data object
    data1 = ConstructionBase.setproperties(data, (value=readval, indices = I))
    # Pass all of these to the applyrule function
    writeval = applyrule(data1, rule, readval, I1)
    # And write its result/s to the cell in the relevent grid/s
    _writecell!(data1, ruletype, wkeys, writeval, I...)
    # We also return the written value
    return writeval
end
@inline function cell_kernel!(data::RuleData, ::Val{<:SetRule}, rule, rkeys, wkeys, I...)
    I1 = _strip_replicates(data, I)
    m = mask(data)
    if !isnothing(m)
        m[I1...] || return nothing
    end
    readval = _readcell(data, rkeys, I...)
    data1 = ConstructionBase.setproperties(data, (value=readval, indices = I))
    # Rules will manually write to grids in `applyrule!`
    applyrule!(data1, rule, readval, I1)
    # In a SetRule there is no return value
    return nothing
end

_strip_replicates(data::RuleData, I) = _strip_replicates(replicates(data), I)
_strip_replicates(::Nothing, I::NTuple) = I
_strip_replicates(::Integer, I::NTuple{N}) where N = ntuple(i -> I[i], Val{N-1}())

# stencil_kernel!
# Runs a rule for the current cell/stencil, when there is no
# row-based optimisation
@inline function stencil_kernel!(
    data::RuleData, hoodgrid::GridData, ruletype::Val{<:NeighborhoodRule}, rule, rkeys, wkeys, I...
)
    rule1 = Stencils.rebuild(rule, unsafe_neighbors(stencil(rule), hoodgrid, CartesianIndex(I)))
    cell_kernel!(data, ruletype, rule1, rkeys, wkeys, I...)
end

# row_kernel!
# Run a NeighborhoodRule rule row by row. When we move along a row by one cell, we 
# access only a single new column of data with the height of 4R, and move the existing
# data in the stencil windows array across by one column. This saves on reads
# from the main array.
function row_kernel!(
    simdata::AbstractSimData, grid::GridData{<:GridMode,<:Tuple{I,J},R}, proc, opt::NoOpt,
    ruletype::Val, rule::Rule, rkeys, wkeys, bi
) where {I,J,R}
    B = 2R
    i = _blocktoind(bi, B)
    i > I && return nothing
    # Loop along the block ROW.
    blocklen = min(I, i + B - 1) - i + 1
    for j = 1:J
        # Loop over the COLUMN of windows covering the block
        for b in 1:blocklen
            stencil_kernel!(simdata, grid, ruletype, rule, rkeys, wkeys, i + b - 1, j)
        end
    end
    return nothing
end


#### Utils

# Convert any GridData to GridData{<:ReadMode}
_to_readonly(data::Tuple) = map(_to_readonly, data)
_to_readonly(data::GridData) = GridData{ReadMode}(data)

# _maybemask!
# mask the source grid with the `mask` array, if it exists
_maybemask!(grids::Union{Tuple,NamedTuple}) = map(_maybemask!, grids)
_maybemask!(grid::GridData) = _maybemask!(grid, proc(grid), mask(grid))
_maybemask!(grid::GridData, proc, mask::Nothing) = nothing
function _maybemask!(grid::GridData, proc::CPU, mask::AbstractArray)
    mv = maskval(grid)
    # `mask` is also a padded StencilArray
    # so we mask the whole thing to take care of the edges
    if isnothing(mv) || iszero(mv)
        source(grid) .*= parent(mask)
    else
        source(grid) .= ((a, m) -> m ? a : mv).(source(grid), parent(mask))
    end
    return grid
end

# _cleardest!
# Clear the desination grid. Only needed with some optimisations.
_cleardest!(grid) = _cleardest!(grid, opt(grid))
_cleardest!(grid, opt) = nothing
