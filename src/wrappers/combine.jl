"""
    Combine{R,W} <: MultiRuleWrapper{R,W}

    Combine(f, rules::Tuple)
    Combine(f, rules...)

Combine the results of multiple independent rules with the function `f`.

`f` combines all values for each grid separately. It must accept a `Tuple`
of length `1` or more

Only `ReturnRule` like `CellRule` and `NeighborhoodRule` can be used in `Combine`.
"""
struct Combine{R,W,F,RL} <: MultiRuleWrapper{R,W}
    f::F
    rules::RL
end
function Combine(f, rules::Tuple{Vararg{ReturnRule}})
    wkeys = Tuple{union(map(k -> _asiterable(_writekeys(k)), rules)...)...} 
    rkeys = Tuple{union(map(k -> _asiterable(_readkeys(k)), rules)...)...}
    Combine{rkeys,wkeys,typeof(f),typeof(rules)}(f, rules)
end
Combine(f, rules::ReturnRule...) = Combine(f, rules)
Combine(rules::ReturnRule...) = Combine(identity, rules)
Combine(rules::Tuple{Vararg{ReturnRule}}) = Combine(identity, rules)

function applyrule(data::AbstractSimData, combinedrule::Combine{R,W}, state, I) where {R,W}
    writestates = map(rules(combinedrule)) do rule
        readstate = _filter_readstate(rule, state)
        applyrule(data, rule, readstate, I)
    end
    return _combinewritestates(combinedrule, writestates)
end

# Combine the write states of multiple rules
@generated function _combinewritestates(combinedrule::Combine{R,W,F,RLS}, writestates::Tuple) where {R,W,F,RLS}
    combineexprs = Expr(:tuple)
    allkeys = map(_writekeys, Tuple(RLS.parameters))
    for key in _writekeys(combinedrule)
        exprargs = Expr(:tuple)
        for (i, rulewritekeys) in enumerate(allkeys)
            if rulewritekeys isa Symbol
                rulewritekeys == key && push!(exprargs.args, :(writestates[$i]))
            else
                j = findfirst(==(key), rulewritekeys)
                isnothing(j) || push!(exprargs.args, :(writestates[$i][$j]))
            end
        end
        push!(combineexprs.args, Expr(:call, :(combinedrule.f), exprargs))
    end
    :($combineexprs)
end
