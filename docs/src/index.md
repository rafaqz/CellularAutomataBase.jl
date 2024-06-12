# DynamicGrids

```@docs
DynamicGrids
```

## Running simulations

```@docs
sim!
resume!
step!
```

## Rulesets

```@docs
AbstractRuleset
Ruleset
```

## Options/Flags

### Boundary conditions

```@docs
BoundaryCondition
Wrap
Remove
Ignore
Reflect
```

### Hardware selection

```@docs
DynamicGrids.Processor
DynamicGrids.CPU
SingleCPU
ThreadedCPU
DynamicGrids.GPU
CuGPU
CPUGPU
```

### Performance optimisation

```@docs
PerformanceOpt
NoOpt
SparseOpt
```

## Rules

```@docs
Rule
DynamicGrids.SetRule
```

### CellRule

```@docs
CellRule
Cell
CopyTo
```

### NeighborhoodRule

```@docs
NeighborhoodRule
Neighbors
Convolution
Life
```

### SetCellRule

```@docs
SetCellRule
SetCell
```

### SetNeighborhoodRule

```@docs
SetNeighborhoodRule
SetNeighbors
```

### SetGridRule

```@docs
SetGridRule
SetGrid
```

### Rule wrappers

```@docs
RuleWrapper
Chain
RunIf
RunAt
```

### Parameter sources

```@docs
ParameterSource
Aux
Grid
DynamicGrids.AbstractDelay
Delay
Frame
Lag
```

### Custom Rule interface and helpers

```@docs
DynamicGrids.applyrule
DynamicGrids.applyrule!
DynamicGrids.modifyrule
isinferred
```

### Methods and objects for use in `applyrule` and/or `modifyrule`

```@docs
get
DynamicGrids.isinbounds
DynamicGrids.inbounds
DynamicGrids.ismasked
DynamicGrids.init
DynamicGrids.aux
DynamicGrids.mask
DynamicGrids.tspan
DynamicGrids.timestep
DynamicGrids.currenttimestep
DynamicGrids.currenttime
DynamicGrids.currentframe
DynamicGrids.AbstractSimData
DynamicGrids.SimData
DynamicGrids.RuleData
DynamicGrids.GridData
DynamicGrids.ReadableGridData
DynamicGrids.WritableGridData
DynamicGrids.AbstractSimSettings
DynamicGrids.SimSettings
```

## Neighborhoods

```@docs
Neighborhood
Moore
VonNeumann
Window
DynamicGrids.AbstractPositionalNeighborhood
Positional
LayeredPositional
```

### Methods for use with neighborhood rules and neighborhoods

```@docs
neighborhood
radius
distances
```

Useful with [`NeighborhoodRule`](@ref):

```@docs
neighbors
```

Useful with [`SetNeighborhoodRule`](@ref):

```@docs
positions
offsets
```

### Convolution kernel neighborhoods

```@docs
AbstractKernelNeighborhood
Kernel
kernel
kernelproduct
```

### Low level use of neighborhoods

```@docs
DynamicGrids.Neighborhoods.readwindow
DynamicGrids.Neighborhoods.unsafe_readwindow
DynamicGrids.Neighborhoods.updatewindow
DynamicGrids.Neighborhoods.unsafe_updatewindow
DynamicGrids.Neighborhoods.pad_axes
DynamicGrids.Neighborhoods.unpad_axes
```

### Generic neighborhood applicators 

These can be used without the full simulation mechanisms, like `broadcast`.

```@docs
DynamicGrids.Neighborhoods.broadcast_neighborhood
DynamicGrids.Neighborhoods.broadcast_neighborhood!
```


## Atomic methods for SetCellRule and SetNeighborhoodRule

Using these methods to modify grid values ensures cell independence,
and also prevent race conditions with [`ThreadedCPU`](@ref) or [`CuGPU`].

```@docs
add!
sub!
min!
max!
and!
or!
xor!
```

## Output

### Output Types and Constructors

```@docs
Output
ArrayOutput
ResultOutput
TransformedOutput
GraphicOutput
REPLOutput
ImageOutput
GifOutput
```

### Renderers

```@docs
Renderer
DynamicGrids.SingleGridRenderer
Image
DynamicGrids.MultiGridRenderer
Layout
SparseOptInspector
```

### Color schemes

Schemes from Colorschemes.jl can be used for the `scheme` argument to `ImageOutput`, 
`Renderer`s. `Greyscale` control over the band of grey used, and is very fast. 
`ObjectScheme` is the default.

```@docs
ObjectScheme
Greyscale
```

### Text labels

```@docs
TextConfig
DynamicGrids.set_default_font
```


### Saving gifs

```@docs
savegif
```

### `Output` interface

These are used for defining your own outputs and `GridProcessors`,
not for general scripting.

```@docs
DynamicGrids.AbstractExtent
DynamicGrids.Extent
DynamicGrids.extent
DynamicGrids.isasync
DynamicGrids.storeframe!
DynamicGrids.isrunning
DynamicGrids.isshowable
DynamicGrids.isstored
DynamicGrids.initialise!
DynamicGrids.finalise!
DynamicGrids.frameindex
```

### `GraphicOutput` interface

Also includes `Output` interface.

```@docs
DynamicGrids.GraphicConfig
DynamicGrids.graphicconfig
DynamicGrids.fps
DynamicGrids.setfps!
DynamicGrids.showframe
DynamicGrids.initialisegraphics
DynamicGrids.finalisegraphics
```

### `ImageOutput` components and interface

Also uses `Output` and `GraphicOutput` interfaces.

```@docs
DynamicGrids.ImageConfig
DynamicGrids.imageconfig
DynamicGrids.showimage
DynamicGrids.render!
DynamicGrids.to_rgb
```

## Custom grid element types

It is common to use `Bool`, `Int` or `Float64` as the contents of a grid.
But a range of object types can be used if they meet the interface criteria.

Immutable, `isbits` objects are usually better and the only type officially to
work - as they are loaded directly in the simulation. Mutable objects,
especially containing pointers, may lead to incorrect stored results, and wont
work at all on GPUs.

For a custom grid element to work, it must have a number of methods defined.

Methods to define are: 

Minimum:

- `zero`: define zero of the object type

Context dependent, and visualisation:
- `oneunit`: define one of the object type
- `isless`: define comparison between two of the objects
- `*`: multiplication by a `Real` scalar.
- `/`: division by a `Real` scalar.
- `+`: addition to another object of the same type
- `-`: subtraction from another object of the same type
- `to_rgb`: optional: return an `ARGB32` to visualise the object as a pixel 

In this example we define a struct with two fields. You will need to determine the
correct behaviours for your own types, but hopefully this will get you started.

```julia
struct MYStruct{A,B}
    a::A
    b::B
end

Base.isless(a::MyStruct, b::MyStruct) = isless(a.a, b.a)
Base.zero(::Type{<:MyStruct{T1,T2}}) where {T1,T2} = MyStruct(zero(T1), zero(T2))
Base.oneunit(::Type{<:MyStruct{T1,T2}}) where {T1,T2} = MyStruct(one(T1), one(T2))

Base.:*(x::MyStruct, x::Number) = MyStruct(x.a * x, x.b * x)
Base.:*(x::Number, x::MyStruct) = MyStruct(x * x.a, x * x.b)
Base.:/(x::MyStruct, x::Number) = MyStruct(x.a / x, x.b / x)
Base.:+(x1::MyStruct, x2::MyStruct) = MyStruct(x1.a + x2.a, x1.b + x2.b)
Base.:-(x1::MyStruct, x2::MyStruct) = MyStruct(x1.a - x2.a, x1.b - x2.b)
```

To generate rgb colors for an `ImageOuput`, you must define `to_rgb`, 
at least for the default `ObjectScheme`, but this can also be done for other 
schemes such as ColorSchemes.jl, or `GreyScale`, by calling `get` on the scheme and a 
`Real` value. Note that the objects will be normalised to values between zero and one
by `minval` and `maxval` scalars prior to this, using the division operators defined 
above. It is preferable to use `minval` and `maxval` over normalising in `to_rgb` - 
as this will not be as flexible for scripting.

```julia
DynamicGrids.to_rgb(::ObjectScheme, obj::MyStruct) = ARGB32(obj.a, obj.b, 0)
DynamicGrids.to_rgb(scheme, obj::MyStruct) = get(scheme, obj.a)
```

See the `test/objectgrids.jl` tests for more details on using complex objects in grids.
