
"""
    ismasked(data, I...)

Check if a cell is masked, using the `mask` array.

Used used internally during simulations to skip masked cells.

If `mask` was not passed to the `Output` constructor or `sim!`
it defaults to `nothing` and `false` is always returned.
"""
ismasked(data::AbstractSimData, I...) = ismasked(mask(data), I...)
ismasked(data::GridData, I...) = ismasked(mask(data), I...)
ismasked(mask::Nothing, I...) = false
ismasked(mask::AbstractArray, I...) = @inbounds !(mask[I...])

abstract type DemoMode end
struct Infer <: DemoMode end
struct Descend <: DemoMode end

"""
    isinferred(output::Output, ruleset::Ruleset)
    isinferred(output::Output, rules::Rule...)

Test if a custom rule is inferred and the return type is correct when
`applyrule` or `applyrule!` is run.

Type-stability can give orders of magnitude improvements in performance.
"""
isinferred(args...) = demo(Infer(), args...)

descendable(args...) = demo(Descend(), args...)

demo(m::DemoMode, output::Output, rules::Rule...) = demo(m, output, rules)
demo(m::DemoMode, output::Output, rules::Tuple) = demo(m, output, StaticRuleset(rules...))
demo(m::DemoMode, output::Output, ruleset::Ruleset) = demo(m, output, StaticRuleset(ruleset))
function demo(m, output::Output, ruleset::StaticRuleset)
    simdata = _updaterules(rules(ruleset), SimData(output, ruleset))
    return _demo(m, simdata)
end
@noinline function demo(m::DemoMode, simdata::AbstractSimData)
    # sd = _updaterules(rules(simdata), simdata)
    _demo(m, simdata)
end
@noinline function _demo(m::DemoMode, simdata::AbstractSimData)
    map(rules(simdata)) do rule
        if rule isa SetGridRule
            return true
        else
            _demo(m, RuleData(simdata), rule, Val{ruletype(rule)}())
        end
    end |> all
end
@noinline function _demo(m::DemoMode, ruledata::RuleData{<:Any,N}, 
   rule, ruletype::Val{<:Union{<:NeighborhoodRule,<:Chain{<:Any,<:Any,<:Tuple{<:NeighborhoodRule,Vararg}}}}
) where N
    hoodgrid = ruledata[stencilkey(rule)]
    rkeys, rgrids = _getreadgrids(rule, ruledata)
    wkeys, wgrids = _getwritegrids(WriteMode, rule, ruledata)
    ruledata = @set ruledata.grids = _combinegrids(rkeys, rgrids, wkeys, wgrids)
    I = ntuple(_ -> 1, Val{N}())
    example_writeval = _example_writeval(wgrids, I)
    if m isa Infer
        writeval = @inferred stencil_kernel!(ruledata, hoodgrid, ruletype, rule, rkeys, wkeys, I...)
        typeof(writeval) == typeof(example_writeval) ||
            error("return type `$(typeof(writeval))` doesn't match grid types `$(typeof(example_writeval))`")
    else
        stencil_kernel!(ruledata, hoodgrid, ruletype, rule, rkeys, wkeys, I...)
    end
    return true
end
@noinline function _demo(m::DemoMode, ruledata::RuleData{<:Any,N}, rule, ruletype::Val{<:CellRule}) where N
    rkeys, rgrids = _getreadgrids(rule, ruledata)
    wkeys, wgrids = _getwritegrids(WriteMode, rule, ruledata)
    ruledata = @set ruledata.grids = _combinegrids(rkeys, rgrids, wkeys, wgrids)
    I = ntuple(_ -> 1, Val{N}())
    example_writeval = _example_writeval(wgrids, I)
    if m isa Infer
        writeval = @inferred cell_kernel!(ruledata, ruletype, rule, rkeys, wkeys, I...)
        typeof(writeval) == typeof(example_writeval) ||
            error("return type `$(typeof(writeval))` doesn't match grid types `$(typeof(example_writeval))`")
    else
        cell_kernel!(ruledata, ruletype, rule, rkeys, wkeys, I...)
    end
    return true
end
@noinline function _demo(m::DemoMode, ruledata::RuleData{<:Any,N}, rule, ruletype::Val{<:SetRule}) where N
    rkeys, rgrids = _getreadgrids(rule, ruledata)
    wkeys, wgrids = _getwritegrids(WriteMode, rule, ruledata)
    ruledata = @set ruledata.grids = _combinegrids(rkeys, rgrids, wkeys, wgrids)
    I = ntuple(_ -> 1, Val{N}())
    if m isa Infer
        @inferred cell_kernel!(ruledata, ruletype, rule, rkeys, wkeys, I...)
    else
        cell_kernel!(ruledata, ruletype, rule, rkeys, wkeys, I...)
    end
    return true
end

_example_writeval(grids::Tuple, I) = Tuple(map(g -> _example_writeval(g, I), grids))
_example_writeval(grid::GridData{<:WriteMode}, I) = grid[I...]

# _zerogrids
# Generate a Vector of zero valued grids
_zerogrids(initgrid::AbstractArray, length) = [zero(initgrid) for f in 1:length]
_zerogrids(initgrids::NamedTuple, length) =
    [map(grid -> zero(grid), initgrids) for f in 1:length]

# _asiterable
# Return some iterable value from a 
# Symbol, Tuple or tuple type
@inline _asiterable(x) = (x,)
@inline _asiterable(x::Symbol) = (x,)
@inline _asiterable(x::Type{<:Tuple}) = x.parameters
@inline _asiterable(x::Tuple) = x
@inline _asiterable(x::AbstractArray) = x

# _astuple
# Wrap a value in a tuple if the matching keys are not a tuple
# we cant just dispatch on state, as it may be meant to be a tuple.
@inline _astuple(rule::Rule, state) = _astuple(_readkeys(rule), state)
@inline _astuple(keys::Tuple, state) = state
@inline _astuple(key, state) = (state,)

# _asnamedtuple => NamedTuple
# Returns a NamedTuple given a NamedTuple or an Array.
# the Array will be called _default_.
@inline _asnamedtuple(x::NamedTuple) = x
@inline _asnamedtuple(x) = (_default_=x,)
@inline function _asnamedtuple(e::AbstractExtent) 
    init_nt = _asnamedtuple(init(e))
    e = @set e.init = init_nt
    pv = _samenamedtuple(init_nt, padval(e))
    return @set e.padval = pv
end

# _samenamedtuple => NamedTuple
# Returns a NamedTuple with length and keys matching the `init` 
# NamedTuple, for another NamedTuple, a Tuple, or a scalar.
@inline _samenamedtuple(init::NamedTuple{K}, x::NamedTuple{K}) where K = x
@noinline _samenamedtuple(init::NamedTuple{K}, x::NamedTuple{J}) where {K,J} = 
    error("Keys $K and $J do not match")
@inline _samenamedtuple(init::NamedTuple{K}, x::Tuple) where K = NamedTuple{K}(x)
@inline _samenamedtuple(init::NamedTuple, x) = map(_ -> x, init) 


# Unwrap a Val or Val type to its internal value
_unwrap(x) = x
_unwrap(::Val{X}) where X = X
_unwrap(::Type{<:Val{X}}) where X = X

@inline _firstgrid(simdata, ::Val{K}) where K = simdata[K]
@inline _firstgrid(simdata, ::Tuple{Val{K},Vararg}) where K = simdata[K]

_replicate_init(init, ::Nothing) = deepcopy(init)
function _replicate_init(init::NamedTuple, replicates::Integer)
    map(init) do layer
        _replicate_init(layer, replicates)
    end
end
function _replicate_init(init::AbstractArray, replicates::Integer)
    A = similar(init, size(init)..., replicates)
    _replicate_init!(A, init, replicates)
    return A
end
function _replicate_init(init::AbstractDimArray, replicates::Integer)
    rep_data = _replicate_init(parent(init), replicates)
    rep_dims = (dims(init)..., Dim{:__replicates__}(NoLookup(Base.OneTo(replicates))))
    return rebuild(init; data=rep_data, dims=rep_dims)
end

function _replicate_init!(A::AbstractArray, init::AbstractArray, replicates::Integer)
    for r in 1:replicates
        A[:, :, r] .= init
    end
end
