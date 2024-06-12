module DynamicGrids
# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end DynamicGrids


import Adapt,
      ConstructionBase,
      FreeTypeAbstraction,
      Reexport,
      REPL

using Colors,
      Crayons,
      DimensionalData,
      FileIO,
      LinearAlgebra,
      KernelAbstractions,
      Stencils,
      Setfield,
      StaticArrays,
      Test,
      UnicodeGraphics

Reexport.@reexport using ModelParameters

import ModelParameters.Flatten

using DimensionalData.LookupArrays, DimensionalData.Dimensions

const DG = DynamicGrids
const DD = DimensionalData

using Base: tail, @propagate_inbounds

import Base: show, getindex, setindex!, lastindex, size, length, push!, append!,
             broadcast, broadcast!, similar, eltype, iterate


export sim!, resume!, step!, savegif, isinferred

export rules

export inbounds, isinbounds 

export gridsize, currentframe, currenttime, currenttimestep, timestep

export add!, sub!, min!, max!, and!, or!, xor!

export Extent, StaticExtent

export Rule, NeighborhoodRule, CellRule, SetCellRule, SetNeighborhoodRule, SetGridRule

export Cell, Neighbors, SetNeighbors, SetCell, Convolution, SetGrid, Life, CopyTo

export RuleWrapper, Chain, RunIf, RunAt, Combine

export AbstractRuleset, Ruleset, StaticRuleset

export AbstractSimData

export Processor, SingleCPU, ThreadedCPU, CPUGPU, CuGPU

export PerformanceOpt, NoOpt, SparseOpt

export BoundaryCondition, Remove, Wrap, Use, Reflect

export ParameterSource, Aux, Grid, Delay, Lag, Frame

export Output, ArrayOutput, ResultOutput, TransformedOutput

export GraphicOutput, REPLOutput, MakieOutput

export ImageOutput, GifOutput

export Renderer, Image, Layout, SparseOptInspector

export TextConfig

export ObjectScheme, Greyscale, Grayscale

export CharStyle, Block, Braile

# From Stencils module
export Stencil, Window, Kernel, Moore, VonNeumann, Positional, Layered

export neighbors, stencil, kernel, kernelproduct, offsets, indices, radius, distances, distance_zones

const DEFAULT_KEY = :_default_

const EXPERIMENTAL = """
        WARNING: This feature is experimental. It may change in future versions, 
        and may not be 100% reliable in all cases. Please file github issues if problems occur.
        """

import Stencils: neighbors, unsafe_neighbors, stencil,
    kernel, kernelproduct, offsets, indices, radius, distances, distance_zones,
    stencil, unsafe_stencil, boundary, padding, source, dest, switch, padval, update_boundary!, add_halo

import Stencils: BoundaryCondition, Padding

include("interface.jl")
include("flags.jl")
include("rules.jl")
include("settings.jl")
include("rulesets.jl")
include("extent.jl")
include("grid.jl")
include("simulationdata.jl")
include("parametersources.jl")
include("gpu.jl")
include("atomic.jl")
include("wrappers/wrappers.jl")
include("wrappers/chain.jl")
include("wrappers/condition.jl")
include("wrappers/combine.jl")
include("outputs/interface.jl")
include("outputs/output.jl")
include("outputs/graphic.jl")
include("outputs/schemes.jl")
include("outputs/textconfig.jl")
include("outputs/render.jl")
include("outputs/image.jl")
include("outputs/makie.jl")
include("outputs/array.jl")
include("outputs/transformed.jl")
include("outputs/repl.jl")
include("outputs/gif.jl")
include("framework.jl")
include("modifyrule.jl")
include("sequencerules.jl")
include("generated.jl")
include("maprules.jl")
include("boundaries.jl")
include("sparseopt.jl")
include("utils.jl")
include("copyto.jl")
include("life.jl")
include("adapt.jl")
include("show.jl")

function __init__()
    global terminal
    terminal = REPL.Terminals.TTYTerminal(get(ENV, "TERM", Base.Sys.iswindows() ? "" : "dumb"), stdin, stdout, stderr)
end

end
