
"""
    sim!(output, rules::Rule...; kw...)
    sim!(output, rules::Tuple{<:Rule,Vararg}; kw...)
    sim!(output, [ruleset::Ruleset=ruleset(output)]; kw...)

Runs the simulation rules over the `output` `tspan`,
writing the destination array to `output` for each time-step.

# Arguments

- `output`: An [`Output`](@ref) to store grids or display them on the screen.
- `ruleset`: A [`Ruleset`](@ref) containing one or more [`Rule`](@ref)s. If the output
  has a `Ruleset` attached, it will be used.

# Keywords

Theses are the taken from the `output` argument by default:

- `init`: optional array or NamedTuple of arrays.
- `mask`: a `Bool` array matching the init array size. `false` cells do not run.
- `aux`: a `NamedTuple` of auxilary data to be used by rules.
- `tspan`: a tuple holding the start and end of the timespan the simulaiton will run for.
- `fps`: the frames per second to display. Will be taken from the output if not passed in.

Theses are the taken from the `ruleset` argument by default:

- `proc`: a [`Processor`](@ref) to specificy the hardware to run simulations on, 
    like [`SingleCPU`](@ref), [`ThreadedCPU`](@ref) or [`CuGPU`](@ref) when 
    KernelAbstractions.jl and a CUDA gpu is available. 
- `opt`: a [`PerformanceOpt`](@ref) to specificy optimisations like
    [`SparseOpt`](@ref) or [`NoOpt`](@ref). Defaults to `NoOpt()`.
- `boundary`: what to do with boundary of grid edges.  Options are [`Remove`](@ref), [`Wrap`](@ref) or [`Reflect`](@ref), defaulting to `Remove()`.
- `cellsize`: the size of cells, which may be accessed by rules.
- `timestep`: fixed timestep where this is required for some rules.
    eg. `Month(1)` or `1u"s"`.

Other:

- `simdata`: a [`SimData`](@ref) object. Keeping it between simulations can reduce memory
  allocation a little, when that is important.
"""
function sim!(output::Output, ruleset::AbstractRuleset=ruleset(output);
    init=init(output),
    mask=mask(output),
    tspan=tspan(output),
    aux=aux(output),
    fps=fps(output),
    boundary=boundary(ruleset),
    proc=proc(ruleset),
    opt=opt(ruleset),
    cellsize=cellsize(ruleset),
    timestep=timestep(ruleset),
    simdata=nothing, 
    kw...
)
    # isrunning(output) && error("Either a simulation is already running in this output, or an error occurred")
    setrunning!(output, true) || error("Could not start the simulation with this output")

    gridsize(init) == gridsize(DG.init(output)) || throw(ArgumentError("init size does not match output init"))

    # Rebuild Extent to allow kwarg alterations
    extent = Extent(; 
        init=_asnamedtuple(init),
        mask=mask,
        aux=aux,
        padval=_asnamedtuple(padval(output)),
        replicates=replicates(output),
        tspan=tspan
    )
    simruleset = Ruleset(rules(ruleset);
        boundary=boundary, proc=proc, opt=opt, cellsize=cellsize, timestep=timestep,
    )
    # Some rules are only valid for a set time-step size.
    step(simruleset) !== nothing && step(simruleset) != step(tspan) &&
        throw(ArgumentError("tspan step $(step(tspan)) must equal rule step $(step(simruleset))"))
    # Set up output
    settspan!(output, tspan)
    # Create or update the combined data object for the simulation
    simdata = SimData(simdata, output, extent, simruleset)
    # Run validation for the rules - they can check if simdata has what they need
    # So error messages happen early, without pages of scrolling.
    _validaterules(ruleset, simdata)
    init_output_grids!(output, init)
    # Set run speed for GraphicOutputs
    setfps!(output, fps)
    # Run any initialisation the output needs to do
    initialise!(output, simdata)
    # Show the first grid
    showframe(output, simdata)
    # Let the init grid be displayed for as long as a normal grid
    maybesleep(output, 1)
    # Run the simulation over simdata and a unitrange we keep 
    # the original ruleset to allow interactive updates to rules.
    # We pass throught the original ruleset as a handle for e.g. 
    # control sliders to update the rules.
    return runsim!(output, simdata, ruleset, 1:lastindex(tspan); kw...)
end
sim!(output::Output, rules::Tuple; kw...) = sim!(output, rules...; kw...)
sim!(output::Output, rules::Rule...; kw...) = sim!(output, Ruleset(rules...; kw...); kw...)

"""
    resume!(output::GraphicOutput, ruleset::Ruleset=ruleset(output); tstop, kw...)

Restart the simulation from where you stopped last time. For arguments see [`sim!`](@ref).
The keyword arg `tstop` can be used to extend the length of the simulation.

# Arguments

- `output`: An [`Output`](@ref) to store grids or display them on the screen.
- `ruleset`: A [`Ruleset`](@ref) containing one ore more [`Rule`](@ref)s.
    These will each be run in sequence.

# Keywords (optional)

- `tstop`: the new stop time for the simulation. Taken from the output length by default.
- `fps`: the frames per second to display. Taken from the output by default.
- `simdata`: a [`SimData`](@ref) object. Keeping it between simulations can improve performance
    when that is important
"""
function resume!(output::GraphicOutput, ruleset::Ruleset=ruleset(output);
        tstop=last(tspan(output)),
        fps=fps(output),
        simdata=nothing,
        kw...
)
    # Check status and arguments
    isrunning(output) && error("A simulation is already running in this output")

    # Calculate new timespan
    new_tspan = first(tspan(output)):step(tspan(output)):tstop
    frame = stoppedframe(output)
    fspan = frame:lastindex(new_tspan)
    settspan!(output, new_tspan)

    # Use the last frame of the existing simulation as the init frame
    if frame <= length(output)
        init = output[frame]
    else
        init = output[1]
    end

    setfps!(output, fps)
    extent = Extent(; init=_asnamedtuple(init), mask=mask(output), aux=aux(output), tspan=new_tspan)
    simdata = SimData(simdata, output, extent, ruleset)
    initialise!(output, simdata)
    setrunning!(output, true) || error("Could not start the simulation with this output")
    return runsim!(output, simdata, ruleset, fspan; kw...)
end

# Simulation runner. Runs a simulation synchonously or asynchonously
# depending on the return value of `isasync(output)` - which may be a
# fixed trait or a field value depending on the output type.

# This allows interfaces with interactive components to update during the simulations.
function runsim!(output, simdata, ruleset, fspan; kw...)
    if isasync(output)
        @async simloop!(output, simdata, ruleset, fspan; kw...)
    else
        simloop!(output, simdata, ruleset, fspan; kw...)
    end
end

# Loop over the frames in `fspan`, running the ruleset and displaying the output.

# Operations on outputs and rulesets are allways mutable and in-place.

# Operations on [`Rule`](@ref)s and [`SimData`](@ref) objects are in a
# functional style, as they are used in inner loops where immutability improves
# performance.
function simloop!(output::Output, simdata, ruleset, fspan; printframe=false)
    printframe && _printframe(simdata, 1)
    # Generate any initialisation data the rules need
    rule_initialisation = initialiserules(simdata)
    # Set the frame timestamp for fps calculation
    settimestamp!(output, first(fspan))
    # Initialise types etc
    simdata = _updatetime(simdata, 1) |> _proc_setup
    # Loop over the simulation
    for f in fspan[2:end]
        printframe && _printframe(simdata, f)
        # Update the current simulation frame and time
        simdata = _updatetime(simdata, f) 
        # Update any Delay parameters
        drules = _setdelays(rules(ruleset), simdata)
        # Run a timestep
        simdata = _step!(simdata, drules)
        # Save/do something with the the current grid
        storeframe!(output, simdata)
        # Let output UI things happen
        yield()
        # Stick to the FPS
        maybesleep(output, f)
        # Exit gracefully
        if !isrunning(output) || f == last(fspan)
            showframe(output, simdata)
            setstoppedframe!(output, f)
            finalise!(output, simdata)
            break
        end
    end
    setrunning!(output, false)
    return output
end

_printframe(simdata, i) = println(stdout, "frame: $i, time: $(tspan(simdata)[i])")

_step!(sd::AbstractSimData, rules) = _updaterules(rules, sd) |> sequencerules!

"""
    step!(sd::AbstractSimData)

Allows stepping a simulation one frame at a time, for a more manual approach
to simulation that `sim!`. This may be useful if other processes need to be run 
between steps, or the simulation is of variable length. `step!` also removes the use
of `Output`s, meaning storing of grid data must be handled manually, if that is 
required. Of course, an output can also be updated manually, using:

```julia
DynmicGrids.storeframe!(output, simdata)
```

Instead of an `Output`, the internal [`SimData`](@ref) objects are used directly, 
and can be defined using a [`Extent`](@ref) object and a [`Ruleset`](@ref).

# Example

```julia
using DynmicGrids, Plots
ruleset = Ruleset(Life(); proc=ThreadedCPU())
extent = Extent(; init=(a=A, b=B), aux=aux, tspan=tspan)
simdata = SimData(extent, ruleset)

# Run a single step, which returns an updated `SimData` object
simdata = step!(simdata)
# Get a view of the grid without padding
grid = DynmicGrids.gridview(simdata[:a])
heatmap(grid)
```

This example returns a `GridData` object for the `:a` grid, which is `<: AbstractAray`.
"""
step!(sd::AbstractSimData; rules=rules(sd)) = step!(sd::AbstractSimData, rules)

step!(sd::AbstractSimData, r1::Rule, rs::Rule...) = step!(sd::AbstractSimData, (r1, rs...))
function step!(sd::AbstractSimData, rules::Tuple)
    _updatetime(sd, currentframe(sd) + 1) |> _proc_setup |> sd -> _step!(sd, rules)
end

# _proc_setup
# Allows different processors to modify the simdata object
# GPU needs this to convert arrays to CuArray
_proc_setup(simdata::AbstractSimData) = _proc_setup(proc(simdata), simdata)
_proc_setup(proc, simdata) = simdata

# function _checkboxed(
# rnd
