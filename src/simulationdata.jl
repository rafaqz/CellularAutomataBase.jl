
"""
    AbstractSimData

Supertype for simulation data objects. Thes hold [`GridData`](@ref), 
[`SimSettings`](@ref) and other objects needed to run the simulation, 
and potentially required from within rules.

An `AbstractSimData` object is accessable in [`applyrule`](@ref) as the first parameter.

Multiple grids can be indexed into using their key if you need to read
from arbitrary locations:

```julia
funciton applyrule(data::AbstractSimData, rule::SomeRule{Tuple{A,B}},W}, (a, b), I) where {A,B,W}
    grid_a = data[A]
    grid_b = data[B]
    ...
end
```

In single-grid simulations `AbstractSimData` objects can be indexed directly as 
if they are a `Matrix`.

## Methods

- `currentframe(data)`: get the current frame number, an `Int`
- `currenttime(data)`: the current frame time, which `isa eltype(tspan)`
- `aux(data, args...)`: get the `aux` data `NamedTuple`, or `Nothing`.
    adding a `Symbol` or `Val{:symbol}` argument will get a field of aux.
- `tspan(data)`: get the simulation time span, an `AbstractRange`.
- `timestep(data)`: get the simulaiton time step.
- `boundary(data)` : returns the [`BoundaryCondition`](@ref) - `Remove` or `Wrap`.
- `padding(data)` : returns the value to use as grid border padding.

These are also available, but you probably shouldn't use them and their behaviour
is not guaranteed in furture versions. Using them will also mean a rule is useful 
only in specific contexts, which is discouraged.

- `settings(data)`: get the simulaitons [`SimSettings`](@ref) object.
- `extent(data)` : get the simulation [`AbstractExtent`](@ref) object.
- `init(data)` : get the simulation init `AbstractArray`/`NamedTuple`
- `mask(data)` : get the simulation mask `AbstractArray`
- `source(data)` : get the `source` grid that is being read from.
- `dest(data)` : get the `dest` grid that is being written to.
- `radius(data)` : returns the `Int` radius used on the grid,
    which is also the amount of border padding.
"""
abstract type AbstractSimData{S,N,G} end

# Getters
extent(d::AbstractSimData) = d.extent
frames(d::AbstractSimData) = d.frames
grids(d::AbstractSimData) = d.grids
auxframe(d::AbstractSimData) = d.auxframe
currentframe(d::AbstractSimData) = d.currentframe

# Forwarded to the Extent object
gridsize(d::AbstractSimData) = size(d)
padval(d::AbstractSimData) = padval(extent(d))
init(d::AbstractSimData) = init(extent(d))
mask(d::AbstractSimData) = mask(extent(d))
aux(d::AbstractSimData, args...) = aux(extent(d), args...)
auxframe(d::AbstractSimData, key) = auxframe(d)[_unwrap(key)]
tspan(d::AbstractSimData) = tspan(extent(d))
timestep(d::AbstractSimData) = step(tspan(d))
radius(d::AbstractSimData) = max(map(radius, grids(d))...)

# Calculated:
# Get the current time for this frame
currenttime(d::AbstractSimData) = tspan(d)[currentframe(d)]
# Get the actual current timestep, e.g. in seconds instead of variable periods like Month
currenttimestep(d::AbstractSimData) = currenttime(d) + timestep(d) - currenttime(d)

# Base methods forwarded to grids NamedTuple
Base.keys(d::AbstractSimData) = keys(grids(d))
Base.keys(d::Type{<:AbstractSimData{<:Any,<:Any,<:NamedTuple{<:K}}}) where K = K
Base.values(d::AbstractSimData) = values(grids(d))
Base.first(d::AbstractSimData) = first(grids(d))
Base.last(d::AbstractSimData) = last(grids(d))
Base.getindex(d::AbstractSimData, key::Symbol) = getindex(grids(d), key)
Base.ndims(d::AbstractSimData{<:Any,N}) where N = N
Base.size(d::AbstractSimData{S}) where S = Tuple(StaticArrays.Size(S))

# Indexing forwarded to the first grid
@propagate_inbounds Base.setindex!(d::AbstractSimData, x, I...) =
    setindex!(first(grids(d)), x, I...)
@propagate_inbounds Base.getindex(d::AbstractSimData, I...) = getindex(first(grids(d)), I...)

# Uptate timestamp
function _updatetime(simdata::AbstractSimData, f::Integer) 
    @set! simdata.currentframe = f
    @set simdata.auxframe = _calc_auxframe(simdata)
end

"""
    SimData <: AbstractSimData

    SimData(extent::AbstractExtent, ruleset::AbstractRuleset)

Simulation dataset to hold all intermediate arrays, timesteps
and frame numbers for the current frame of the simulation.

Additional methods not found in [`AbstractSimData`](@ref):

- `rules(d::SimData)` : get the simulation rules.
- `ruleset(d::SimData)` : get the simulation [`AbstractRuleset`](@ref).
"""
struct SimData{S<:Tuple,N,G<:NamedTuple{<:Any,<:Tuple{<:GridData,Vararg{<:GridData}}},E,RS,F,CF,AF} <: AbstractSimData{S,N,G}
    grids::G
    extent::E
    ruleset::RS
    frames::F
    currentframe::CF
    auxframe::AF
end
function SimData{S,N}(
    grids::G, extent::E, ruleset::RS, frames::F, currentframe::CF, auxframe::AF
) where {S,N,G,E,RS,F,CF,AF}
    SimData{S,N,G,E,RS,F,CF,AF}(grids, extent, ruleset, frames, currentframe, auxframe)
end
SimData(o, ruleset::AbstractRuleset) = SimData(o, extent(o), ruleset)
SimData(o, r1::Rule, rs::Rule...) = SimData(o, extent(o), Ruleset(r1, rs...))

# When no simdata is passed in, create new SimData
function SimData(::Nothing, output, extent::AbstractExtent, ruleset::AbstractRuleset)
    SimData(output, extent, ruleset)
end
# Initialise a AbstractSimData object with a new `Extent` and `Ruleset`.
function SimData(
    simdata::SimData, output, extent::AbstractExtent, ruleset::AbstractRuleset
)
    (replicates(simdata) == replicates(output) == replicates(extent)) || 
        throw(ArgumentError("`simdata` must have same numver of replicates as `output`"))

    @assert simdata.extent == StaticExtent(extent)
    @set! simdata.ruleset = StaticRuleset(ruleset)
    if hasdelay(rules(ruleset)) 
        isstored(output) || _not_stored_delay_error()
        @set! simdata.frames = frames(output) 
    end
    return simdata
end

function SimData(o, extent::AbstractExtent, ruleset::AbstractRuleset)
    frames_ = if hasdelay(rules(ruleset)) 
        isstored(o) || _notstorederror()
        frames(o) 
    else
        nothing 
    end
    return SimData(extent, ruleset, frames_)
end
# Convenience constructors
SimData(extent::AbstractExtent, r1::Rule, rs::Rule...) = SimData(extent, (r1, rs...))
SimData(extent::AbstractExtent, rs::Tuple{<:Rule,Vararg}) = SimData(extent, Ruleset(rs))
# Convert grids in extent to NamedTuple
function SimData(extent::AbstractExtent, ruleset::AbstractRuleset, frames=nothing)
    nt_extent = _asnamedtuple(extent)
    SimData(nt_extent, ruleset, frames) 
end
function SimData(extent::AbstractExtent{<:NamedTuple}, ruleset::AbstractRuleset, frames=nothing)
    # Calculate the stencil array for each grid
    grids = _buildgrids(extent, ruleset)
    # Construct the SimData for each grid
    SimData(grids, extent, ruleset, frames)
end
function SimData(
    grids::G, extent::AbstractExtent, ruleset::AbstractRuleset, frames
) where {G<:Union{<:NamedTuple{<:Any,<:Tuple{<:GridData,Vararg}},<:GridData}}
    currentframe = 1; auxframe = nothing
    S = Tuple{size(extent)...}
    N = ndims(extent)
    # SimData is isbits-only, so use Static versions
    s_extent = StaticExtent(extent)
    s_ruleset = StaticRuleset(ruleset)
    SimData{S,N}(grids, s_extent, s_ruleset, frames, currentframe, auxframe)
end

# Build the grids for the simulation from the extent, ruleset, init and padding
function _buildgrids(extent::AbstractExtent{<:NamedTuple{Keys}}, ruleset) where Keys
    S = Val{Tuple{size(extent)...}}()
    radii = map(k-> Val{get(radius(ruleset), k, 0)}(), Keys)
    radii = NamedTuple{Keys}(radii)
    _buildgrids(extent, ruleset, S, radii)
end
function _buildgrids(extent, ruleset, s::Val, radii::NamedTuple)
    map(radii, init(extent), padval(extent)) do r, in, pv
        _buildgrids(extent, ruleset, s, r, in, pv)
    end
end
function _buildgrids(extent, ruleset, ::Val{S}, ::Val{R}, init, padval) where {S,R}
    stencil = Window{R}() 
    padding = Halo{:out}() # We always pad out in DynamicGrids - it should pay back for multiple time steps
    bc = _update_padval(boundary(ruleset), padval)
    data = _replicate_init(init, replicates(extent))
    m = if isnothing(mask(extent))
        nothing
    else
        Stencils.StencilArray(mask(extent), stencil; padding, boundary=bc)
    end
    GridData{ReadMode,S,R}(
        data, stencil, bc, padding, proc(ruleset), opt(ruleset), m, padval, replicates(extent)
    )
end

_update_padval(::Wrap, padval) = Wrap()
_update_padval(::Remove, padval) = Remove(padval)
_update_padval(::Use, padval) = Use()

ConstructionBase.constructorof(::Type{<:SimData{S,N}}) where {S,N} = SimData{S,N}

# Getters
ruleset(d::SimData) = d.ruleset
rules(d::SimData) = rules(ruleset(d))
boundary(d::SimData) = boundary(ruleset(d))
proc(d::SimData) = proc(ruleset(d))
opt(d::SimData) = opt(ruleset(d))
settings(d::SimData) = settings(ruleset(d))
replicates(d::SimData) = replicates(extent(d))

"""
    RuleData <: AbstractSimData

    RuleData(extent::AbstractExtent, settings::SimSettings)

[`AbstractSimData`](@ref) object that is passed to rules. 
Basically a trimmed-down version of [`SimData`](@ref).

The simplified object actually passed to rules with the current design.

Passing a smaller object than `SimData` to rules leads to faster GPU compilation.
"""
struct RuleData{S<:Tuple,N,G<:NamedTuple,E,Se,F,CF,AF,R,V,I} <: AbstractSimData{S,N,G}
    grids::G
    extent::E
    settings::Se
    frames::F
    currentframe::CF
    auxframe::AF
    replicates::R
    value::V
    indices::I
end
function RuleData{S,N}(
    grids::G, extent::E, settings::Se, frames::F, currentframe::CF, auxframe::AF, replicates::Re, value::V, indices::I
) where {S,N,G,E,Se,F,CF,AF,Re,V,I}
    RuleData{S,N,G,E,Se,F,CF,AF,Re,V,I}(grids, extent, settings, frames, currentframe, auxframe, replicates, value, indices)
end
function RuleData(d::AbstractSimData{S,N};
    grids=grids(d), 
    extent=extent(d), 
    settings=settings(d),
    frames=frames(d), 
    currentframe=currentframe(d), 
    auxframe=auxframe(d), 
    replicates=replicates(d),
    value=nothing, 
    indices=nothing,
) where {S,N}
    RuleData{S,N}(
        grids, extent, settings, frames, currentframe, auxframe, replicates, value, indices
    )
end
# Thin down the aux data for just this rule.
# This may help keep GPU kernel parameter size under limits
function RuleData(d::AbstractSimData, rule::Rule)
    aux_in_rule = Flatten.flatten(rule, Aux)
    rule_aux_keys = map(aux_in_rule) do aux
        _unwrap(aux)
    end
    if isnothing(aux(d))
        if length(rule_aux_keys) > 0
            throw(ArgumentError("no aux data available: expected aux with keys $rule_aux_keys"))
        end
        return RuleData(d)
    end
    rule_aux = aux(d)[rule_aux_keys]
    # Update the extent to have less aux
    rule_extent = ConstructionBase.setproperties(extent(d), (; aux=rule_aux))
    # Use the thinned rule extent in RuleData
    return RuleData(d; extent=rule_extent)
end

function Base.getindex(d::RuleData, key::Symbol) 
    grid = getindex(grids(d), key)
    @set grid.indices = d.indices
end

ConstructionBase.constructorof(::Type{<:RuleData{S,N}}) where {S,N} = RuleData{S,N}

# Getters
settings(d::RuleData) = d.settings
boundary(d::RuleData) = boundary(settings(d))
proc(d::RuleData) = proc(settings(d))
opt(d::RuleData) = opt(settings(d))
replicates(d::RuleData) = d.replicates
indices(d::RuleData) = d.indices
