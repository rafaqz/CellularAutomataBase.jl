"""
    Chain(rules...)
    Chain(rules::Tuple)

`Chain`s allow chaining rules together to be completed in a single processing step,
without intermediate reads or writes from grids.

They are potentially compiled together into a single function call, especially if you
use `@inline` on all `applyrule` methods. `Chain` can hold either all [`CellRule`](@ref)
or [`NeighborhoodRule`](@ref) followed by [`CellRule`](@ref).

[`SetRule`](@ref) can't be used in `Chain`, as it doesn't have a return value.

![Chain rule diagram](https://raw.githubusercontent.com/cesaraustralia/DynamicGrids.jl/media/Chain.png)
"""
struct Chain{R,W,T<:Tuple{Vararg{<:ReturnRule}}} <: MultiRuleWrapper{R,W}
    rules::T
end
Chain(rules...) = Chain(rules)
function Chain(rules::Tuple)
    rkeys = Tuple{union(map(k -> _asiterable(_readkeys(k)), rules)...)...}
    wkeys = Tuple{union(map(k -> _asiterable(_writekeys(k)), rules)...)...}
    Chain{rkeys,wkeys,typeof(rules)}(rules)
end

@inline function Stencils.rebuild(chain::Chain{R,W}, win) where {R,W}
    rules = (Stencils.rebuild(chain[1], win), tail(chain.rules)...)
    Chain{R,W,typeof(rules)}(rules)
end

@generated function applyrule(data::RuleData, chain1::Chain{R,W,T}, state1, index) where {R,W,T}
    expr = Expr(:block)
    nrules = length(T.parameters)
    for i in 1:nrules
        # Variables are numbered to make debugging type stability easier
        state = Symbol("state$i")
        nextstate = Symbol("state$(i+1)") 
        rule = Symbol("rule$i") 
        read = Symbol("read$i") 
        write = Symbol("write$i") 
        chain = Symbol("chain$i") 
        nextchain = Symbol("chain$(i+1)") 
        rule_expr = quote
            $rule = $chain[1]
            # Get the state needed by this rule
            $read = _filter_readstate($rule, $state)
            # Run the rule
            $write = applyrule(data, $rule, $read, index)
            # Create new state with the result and state from other rules
            $nextstate = _update_chainstate($rule, $state, $write)
            $nextchain = tail($chain)
        end
        push!(expr.args, rule_expr)
    end
    laststate = Symbol("state$(nrules+1)")
    push!(expr.args, :(_filter_writestate(data, chain1, $laststate)))
    expr
end

# Merge new state with previous state.
# Returns a new `NamedTuple` with all keys having the most recent state
@generated function _update_chainstate(rule::Rule{R,W}, state::NamedTuple{K,V}, writestate) where {R,W,K,V}
    expr = Expr(:tuple)
    writekeys = W isa Symbol ? (W,) : W.parameters
    keys = (union(K, writekeys)...,)
    for (i, k) in enumerate(keys)
        if k in writekeys
            for (j, wkey) in enumerate(writekeys)
                if k == wkey
                    if W isa Symbol
                        push!(expr.args, :(writestate))
                    else
                        push!(expr.args, :(writestate[$j]))
                    end
                end
            end
        else
            push!(expr.args, :(state[$i]))
        end
    end
    :(NamedTuple{$keys}($expr))
end
