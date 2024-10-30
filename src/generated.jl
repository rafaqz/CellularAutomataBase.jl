# _vals2syms => Union{Symbol,Tuple}
# Must be at the top of the file for world age problems
# Get symbols from a Val or Tuple type
@inline _vals2syms(x::Type{<:Tuple}) = map(v -> _vals2syms(v), x.parameters)
@inline _vals2syms(::Type{<:Val{X}}) where X = X

# Low-level generated functions for working with grids

# _initialise_windows => NTuple{N,SMatrix}
# Generate an SArray from the main array
@generated function _initialise_windows(src::AbstractArray{T}, ::Val{R}, i, j) where {T,R}
    B = 2R; S = 2R + 1; L = S^2
    columns = []
    # This column is never used, so fill with zeros 
    zerocol = Expr[]
    for r in 1:2B
        push!(zerocol, :(zero(T)))
    end
    push!(columns, zerocol)
    for c in 0:S-2
        newcol = Expr[]
        for r in 0:2B-1
            push!(newcol, :(@inbounds src[i + $r, j + $c]))
        end
        push!(columns, newcol)
    end
    newwindows = Expr(:tuple)
    for b in 1:B
        winvals = Expr(:tuple)
        for c in 1:S, r in b:b+B
            exp = columns[c][r]
            push!(winvals.args, exp)
        end
        push!(newwindows.args, :(SArray{Tuple{$S,$S},$T,2,$L}($winvals)))
    end
    return quote
        return $newwindows
    end
end

# _slide_windows => NTuple{N,SMatrix}
# Generate a tuple of SArrays from the main array and the previous SArrays
@generated function _slide_windows(
    windows::Tuple, src::AbstractArray{T}, ::Val{R}, i, j
) where {T,R}
    B = 2R; S = 2R + 1; L = S^2
    newvals = Expr[]
    for n in 0:2B-1
        push!(newvals, :(@inbounds src[i + $n, j + 2R]))
    end
    newwindows = Expr(:tuple)
    for b in 1:B
        winvals = Expr(:tuple)
        for n in S+1:L
            push!(winvals.args, :(@inbounds windows[$b][$n]))
        end
        for n in b:b+B
            push!(winvals.args, newvals[n])
        end
        push!(newwindows.args, :(SArray{Tuple{$S,$S},$T,2,$L}($winvals)))
    end
    return quote
        return $newwindows
    end
end

# _readcell
# Returns a single value or NamedTuple of values.
# This occurs for every cell for every rule, so has to be very fast.
@generated function _readcell(data, ::K, I...) where K<:Tuple
    expr = Expr(:tuple)
    keys = map(_unwrap, Tuple(K.parameters))
    for (i, k) in enumerate(keys)
        push!(expr.args, :(@inbounds source(data[$(QuoteNode(k))])[add_halo(data[$(QuoteNode(k))], I)...]))
    end
    return quote
        keys = $keys
        vals = $expr
        NamedTuple{keys,typeof(vals)}(vals)
    end
end
@inline function _readcell(data::AbstractSimData, ::Val{K}, I...) where K
    grid = data[K]
    @inbounds source(grid)[add_halo(grid, I)...]
end

# _writecell! => nothing
# Write values to grid/s at index `I`. 
# This occurs for every cell for every rule, so has to be very fast.
@generated function _writecell!(
    data, ::Val{<:CellRule}, wkeys::K, vals::Union{Tuple,NamedTuple}, I...
) where K<:Tuple
    expr = Expr(:block)
    keys = map(_unwrap, Tuple(K.parameters))
    for (i, k) in enumerate(keys) 
        # MUST write to source(grid) - CellRule doesn't switch grids
        push!(expr.args, :(@inbounds source(data[$(QuoteNode(k))])[add_halo(data[$(QuoteNode(k))], I)...] = vals[$i]))
    end
    push!(expr.args, :(nothing))
    return expr
end
@inline function _writecell!(data, ::Val{<:CellRule}, wkeys::Val{K}, val, I...) where K
    # MUST write to source(grid) - CellRule doesn't switch grids
    grid = data[K]
    @inbounds source(grid)[add_halo(grid, I)...] = val
    return nothing
end
@generated function _writecell!(
    data, ::Val, wkeys::K, vals::Union{Tuple,NamedTuple}, I...
) where K<:Tuple
    expr = Expr(:block)
    keys = map(_unwrap, Tuple(K.parameters))
    for (i, k) in enumerate(keys) 
        # MUST write to dest(grid) here, not grid K
        # setindex! has overrides for the grid
        push!(expr.args, :(@inbounds dest(data[$(QuoteNode(k))])[add_halo(data[$(QuoteNode(k))], I)...] = vals[$i]))
    end
    push!(expr.args, :(nothing))
    return expr
end
@inline function _writecell!(data, ::Val, wkeys::Val{K}, val, I...) where K
    # MUST write to dest(grid) here, not grid K
    # setindex! has overrides for the grid
    grid = data[K]
    @inbounds dest(grid)[add_halo(grid, I)...] = val
    return nothing
end

# _getgrids => Union{GridData{<:ReadMode},Tuple{GridData{<:ReadMode},Vararg}}
# Retrieves `GridData` from a `SimData` object to match the requirements of a `Rule`.
# Returns a `Tuple` holding the key or `Tuple` of keys, and grid or `Tuple` of grids.
@generated function _getreadgrids(::Rule{R,W}, data::AbstractSimData) where {R<:Tuple,W}
    Expr(:tuple,
        Expr(:tuple, (:(Val{$(QuoteNode(key))}()) for key in R.parameters)...),
        Expr(:tuple, (:(GridData{ReadMode}(data[$(QuoteNode(key))])) for key in R.parameters)...),
    )
end
@generated function _getreadgrids(::Rule{R,W}, data::AbstractSimData) where {R,W}
    :((Val{$(QuoteNode(R))}(), GridData{ReadMode}(data[$(QuoteNode(R))])))
end

# _getwritegrids => Union{GridData{<:WriteMode},Tuple{GridData{<:WriteMode},Vararg}}
# Retrieves `GridData` from a `SimData` object to match the requirements of a `Rule`.
# Returns a `Tuple` holding the key or `Tuple` of keys, and grid or `Tuple` of grids.
@generated function _getwritegrids(
    ::Type{Mode}, ::Rule{R,W}, data::AbstractSimData
) where {Mode<:WriteMode,R,W<:Tuple}
    Expr(:tuple,
        Expr(:tuple, (:(Val{$(QuoteNode(key))}()) for key in W.parameters)...),
        Expr(:tuple, (:(GridData{Mode}(data[$(QuoteNode(key))])) for key in W.parameters)...),
    )
end
@generated function _getwritegrids(
    ::Type{Mode}, ::Rule{R,W}, data::AbstractSimData
) where {Mode<:WriteMode,R,W}
    :((Val{$(QuoteNode(W))}(), GridData{Mode}(data[$(QuoteNode(W))])))
end


# _combinegrids => AbstractSimData
# Combine grids into a new NamedTuple of grids depending
# on the read and write keys required by a rule.
@inline function _combinegrids(data::AbstractSimData, wkeys, wgrids)
    allkeys = map(Val, keys(data))
    allgrids = values(data)
    _combinegrids(data, allkeys, allgrids, wkeys, wgrids)
end
@inline function _combinegrids(data::AbstractSimData, rkeys, rgrids, wkeys, wgrids)
    @set data.grids = _combinegrids(rkeys, rgrids, wkeys, wgrids)
end
@inline function _combinegrids(rkey, rgrids, wkey, wgrids)
    _combinegrids((rkey,), (rgrids,), (wkey,), (wgrids,))
end
@inline function _combinegrids(rkey, rgrids, wkeys::Tuple, wgrids::Tuple)
    _combinegrids((rkey,), (rgrids,), wkeys, wgrids)
end
@inline function _combinegrids(rkeys::Tuple, rgrids::Tuple, wkey, wgrids)
    _combinegrids(rkeys, rgrids, (wkey,), (wgrids,))
end
@generated function _combinegrids(rkeys::Tuple, rgrids::Tuple, wkeys::Tuple, wgrids::Tuple)
    rkeys = _vals2syms(rkeys)
    wkeys = _vals2syms(wkeys)
    keysexp = Expr(:tuple, QuoteNode.(wkeys)...)
    dataexp = Expr(:tuple, :(wgrids...))

    for (i, key) in enumerate(rkeys)
        if !(key in wkeys)
            push!(dataexp.args, :(rgrids[$i]))
            push!(keysexp.args, QuoteNode(key))
        end
    end

    return quote
        keys = $keysexp
        vals = $dataexp
        NamedTuple{keys,typeof(vals)}(vals)
    end
end

# _replacegrids => AbstractSimData
# Replace grids in a NamedTuple with new grids where required.
function _replacegrids(data::AbstractSimData, newkeys, newgrids)
    @set data.grids = _replacegrids(grids(data), newkeys, newgrids)
end
@generated function _replacegrids(allgrids::NamedTuple, newkeys::Tuple, newgrids::Tuple)
    newkeys = map(_unwrap, newkeys.parameters)
    allkeys = allgrids.parameters[1]
    expr = Expr(:tuple)
    for key in allkeys
        if key in newkeys
            i = findfirst(k -> k == key, newkeys)
            push!(expr.args, :(newgrids[$i]))
        else
            push!(expr.args, :(allgrids.$key))
        end
    end

    return quote
        vals = $expr
        NamedTuple{$allkeys,typeof(vals)}(vals)
    end
end
@generated function _replacegrids(allgrids::NamedTuple, newkey::Val, newgrid::GridData)
    newkey = _unwrap(newkey)
    allkeys = allgrids.parameters[1]
    expr = Expr(:tuple)
    for key in allkeys
        if key == newkey
            push!(expr.args, :(newgrid))
        else
            push!(expr.args, :(allgrids.$key))
        end
    end

    return quote
        vals = $expr
        NamedTuple{$allkeys,typeof(vals)}(vals)
    end
end
