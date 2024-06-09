"""
    AbstractRuleset <: ModelParameters.AbstractModel

Abstract supertype for [`Ruleset`](@ref) objects and variants.
""" 
abstract type AbstractRuleset <: AbstractModel end

# Getters
ruleset(rs::AbstractRuleset) = rs
rules(rs::AbstractRuleset) = rs.rules
settings(rs::AbstractRuleset) = rs.settings
boundary(rs::AbstractRuleset) = boundary(settings(rs))
proc(rs::AbstractRuleset) = proc(settings(rs))
opt(rs::AbstractRuleset) = opt(settings(rs))
cellsize(rs::AbstractRuleset) = cellsize(settings(rs))
timestep(rs::AbstractRuleset) = timestep(settings(rs))
function radius(set::AbstractRuleset) 
    if isempty(rules(set))
        NamedTuple()
    else
        radius(rules(set))
    end
end

Base.step(rs::AbstractRuleset) = timestep(rs)

# ModelParameters interface
Base.parent(rs::AbstractRuleset) = rules(rs)
ModelParameters.setparent!(rs::AbstractRuleset, rules) = rs.rules = rules
ModelParameters.setparent(rs::AbstractRuleset, rules) = @set rs.rules = rules

"""
    Rulseset <: AbstractRuleset

    Ruleset(rules...; kw...)
    Ruleset(rules, settings)

A container for holding a sequence of `Rule`s and simulation
details like boundary handing and optimisation.
Rules will be run in the order they are passed, ie. `Ruleset(rule1, rule2, rule3)`.

# Keywords

- `proc`: a [`Processor`](@ref) to specificy the hardware to run simulations on, 
    like [`SingleCPU`](@ref), [`ThreadedCPU`](@ref) or [`CuGPU`](@ref) when 
    KernelAbstractions.jl and a CUDA gpu is available. 
- `opt`: a [`PerformanceOpt`](@ref) to specificy optimisations like
    [`SparseOpt`](@ref). Defaults to [`NoOpt`](@ref).
- `boundary`: what to do with boundary of grid edges.
    Options are `Remove()` or `Wrap()`, defaulting to [`Remove`](@ref).
- `cellsize`: size of cells.
- `timestep`: fixed timestep where this is required for some rules. 
    eg. `Month(1)` or `1u"s"`.
"""
mutable struct Ruleset{S} <: AbstractRuleset
    # Rules in Ruleset are intentionally not type-stable.
    # But they are when rebuilt in a StaticRuleset later
    rules::Tuple{Vararg{Rule}}
    settings::S
end
Ruleset(rule1, rules::Rule...; kw...) = Ruleset((rule1, rules...); kw...)
Ruleset(rules::Tuple; kw...) = Ruleset(rules, SimSettings(; kw...))
Ruleset(rs::AbstractRuleset) = Ruleset(rules(rs), settings(rs))
function Ruleset(; rules=(), settings=nothing, kw...) 
    settings1 = settings isa Nothing ? SimSettings(; kw...) : settings
    Ruleset(rules, settings1)
end

struct StaticRuleset{R<:Tuple,S} <: AbstractRuleset
    rules::R
    settings::S
end
StaticRuleset(rule1, rules::Rule...; kw...) = StaticRuleset((rule1, rules...); kw...)
StaticRuleset(rules::Tuple; kw...) = StaticRuleset(rules, SimSettings(; kw...))
StaticRuleset(rs::AbstractRuleset) = StaticRuleset(rules(rs), settings(rs))
function StaticRuleset(; rules=(), settings=nothing, kw...) 
    settings1 = settings isa Nothing ? SimSettings(; kw...) : settings
    StaticRuleset(rules, settings1)
end

const SRuleset = StaticRuleset
