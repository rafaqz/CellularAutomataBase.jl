const MASKCOL = ARGB32(0.5)
const ZEROCOL = ARGB32(0.3)

"""
    Renderer

Abstract supertype for objects that convert a frame of the simulation into an `ARGB32`
image for display. Frames may be a single grid or a `NamedTuple` of multiple grids.
"""
abstract type Renderer end

imagesize(p::Renderer, init::NamedTuple) = imagesize(p, first(init))
imagesize(::Renderer, init::AbstractArray) = size(init)

_allocimage(p::Renderer, init) = fill(ARGB32(0), imagesize(p, init)...)

function render!(o::ImageOutput, data::AbstractSimData)
    render!(o, data, grids(data))
end
function render!(o::ImageOutput, data::AbstractSimData, grids)
    render!(imagebuffer(o), renderer(o), o, data, grids)
end

"""
    SingleGridRenderer <: Renderer

Abstract supertype for [`Renderer`](@ref)s that convert a single grid 
into an image array.

The first grid will be displayed if a `SingleGridRenderer` is
used with a `NamedTuple` of grids.
"""
abstract type SingleGridRenderer <: Renderer end

function render!(
    imagebuffer, ig::SingleGridRenderer, o::ImageOutput, data::AbstractSimData, 
    grids::NamedTuple;
    name=string(first(keys(grids))), time=currenttime(data)
)
    render!(imagebuffer, ig, o, data, first(grids); name=name, time=time)
end
function render!(
    imagebuffer, ig::SingleGridRenderer, o::ImageOutput, data::AbstractSimData, 
    grids::NamedTuple{(:_default_,)};
    name=nothing, time=currenttime(data)
)
    render!(imagebuffer, ig, o, data, first(grids); name=name, time=time)
end
function render!(
    imagebuffer, ig::SingleGridRenderer, o::ImageOutput, 
    data::AbstractSimData{Y,X}, grid::AbstractArray;
    name=nothing, time=currenttime(data), minval=minval(o), maxval=maxval(o),
) where {Y,X}
    for j in 1:X, i in 1:Y
        @inbounds val = grid[i, j]
        pixel = to_rgb(cell_to_pixel(ig, mask(o), minval, maxval, data, val, (i, j)))
        @inbounds imagebuffer[i, j] = pixel
    end
    _rendertext!(imagebuffer, textconfig(o), name, time)
    return imagebuffer
end

"""
    Image <: SingleGridRenderer

    Image(; scheme=ObjectScheme(), zerocolor=nothing, maskcolor=nothing)

Converts output grids to a colorsheme.

# Keywords

- `scheme`: a ColorSchemes.jl colorscheme, [`ObjectScheme`](@ref) or object that defines
    `Base.get(obj, val)` and returns a `Color` or a value that can be converted to `Color`
    using `ARGB32(val)`.
- `zerocolor`: a `Col` to use when values are zero, or `nothing` to ignore.
- `maskcolor`: a `Color` to use when cells are masked, or `nothing` to ignore.
"""
struct Image{S,Z,M} <: SingleGridRenderer
    scheme::S
    zerocolor::Z
    maskcolor::M
end
Image(scheme, zerocolor=ZEROCOL) = Image(scheme, zerocolor, MASKCOL)
Image(; scheme=ObjectScheme(), zerocolor=ZEROCOL, maskcolor=MASKCOL) =
    Image(scheme, zerocolor, maskcolor)

scheme(p::Image) = p.scheme
zerocolor(p::Image) = p.zerocolor
maskcolor(p::Image) = p.maskcolor

# Show colorscheme in Atom etc
Base.show(io::IO, m::MIME"image/svg+xml", p::Image) = show(io, m, scheme(p))

@inline function cell_to_pixel(ig::Image, mask, minval, maxval, data::AbstractSimData, val, I)
    if !(maskcolor(ig) isa Nothing) && ismasked(mask, I...)
        to_rgb(maskcolor(ig))
    else
        normval = normalise(val, minval, maxval)
        if !(zerocolor(ig) isa Nothing) && normval == zero(typeof(normval))
            to_rgb(zerocolor(ig))
        elseif normval isa Number && isnan(normval)
            zerocolor(ig) isa Nothing ? to_rgb(scheme(ig), 0) : to_rgb(zerocolor(ig))
        else
            to_rgb(scheme(ig), normval)
        end
    end
end

"""
    SparseOptInspector()

A [`Renderer`](@ref) that checks [`SparseOpt`](@ref) visually.
Cells that do not run show in gray. Errors show in red, but if they do there's a bug.
"""
struct SparseOptInspector <: SingleGridRenderer end

function cell_to_pixel(p::SparseOptInspector, mask, minval, maxval, data::AbstractSimData, val, I::Tuple)
    opt(data) isa SparseOpt || error("Can only use SparseOptInspector with SparseOpt grids")
    r = radius(first(grids(data)))
    blocksize = 2r
    blockindex = _indtoblock.((I[1] + r,  I[2] + r), blocksize)
    normedval = normalise(val, minval, maxval)
    status = sourcestatus(first(data))
    # This is done at the start of the next frame, so wont show up in
    # the image properly. So do it preemtively?
    _wrapstatus!(status)
    if status[blockindex...]
        if normedval > 0
            to_rgb(normedval)
        else
            to_rgb((0.0, 0.0, 0.0))
        end
    elseif normedval > 0
        to_rgb((1.0, 0.0, 0.0)) # This (a red cell) would mean there is a bug in SparseOpt
    else
        to_rgb((0.5, 0.5, 0.5))
    end
end

"""
    MultiGridRenderer <: Renderer

Abstract type for `Renderer`s that convert a frame containing multiple 
grids into a single image.
"""
abstract type MultiGridRenderer <: Renderer end

"""
    Layout <: MultiGridRenderer

    Layout(layout::Array, renderer::Matrix)

Layout allows displaying multiple grids in a block layout, by specifying a
layout matrix and a list of [`Image`](@ref)s to be run for each.

# Arguments

- `layout`: A Vector or Matrix containing the keys or numbers of grids in the locations to
    display them. `nothing`, `missing` or `0` values will be skipped.
- `renderer`: tuple of [`Image`](@ref), one for each grid in the simulation.
    Can be `nothing` or any other value for grids not in layout.
"""
Base.@kwdef struct Layout{L<:AbstractMatrix,P} <: MultiGridRenderer
    layout::L
    renderer::P
    Layout(layouts::L, renderer::P) where {L,P} = begin
        imgens = map(_asrenderer, renderer)
        new{L,typeof(imgens)}(layouts, imgens)
    end
end
# Convenience constructor to convert Vector input to a column Matrix
function Layout(layout::AbstractVector, renderer)
    Layout(reshape(layout, length(layout), 1), renderer)
end

_asrenderer(p::Renderer) = p
_asrenderer(x) = Image(x)

layout(p::Layout) = p.layout
renderer(p::Layout) = p.renderer
imagesize(p::Layout, init::NamedTuple) = size(first(init)) .* size(p.layout)

function render!(
    imagebuffer, l::Layout, o::ImageOutput, data::AbstractSimData, grids::NamedTuple
)
    ngrids = length(grids)
    if !(minval(o) isa Nothing)
        length(minval(o)) == ngrids || _wronglengtherror(minval, ngrids, length(minval(o)))
    end
    if !(maxval(o) isa Nothing)
        length(maxval(o)) == ngrids || _wronglengtherror(maxval, ngrids, length(maxval(o)))
    end

    grid_ids = layout(l)
    # Loop over the layout matrix
    for i in 1:size(grid_ids, 1), j in 1:size(grid_ids, 2)
        grid_id = grid_ids[i, j]
        # Accept symbol keys and numbers, skip missing/nothing/0
        (ismissing(grid_id) || grid_id === nothing || grid_id == 0)  && continue
        n = if grid_id isa Symbol
            found = findfirst(k -> k === grid_id, keys(grids))
            found === nothing && _grididnotinkeyserror(grid_id, grids)
            found
        else
            grid_id
        end
        I, J = map((i, j), gridsize(data)) do k, s
            (k - 1) * s + 1:k * s
        end
        # Run image renderers for section
        render!(
            view(imagebuffer, I, J), renderer(l)[n], o, data, grids[n];
            name=string(keys(grids)[n]), time=nothing,
            minval=_valn(minval(o), n), maxval=_valn(maxval(o), n),
        )
    end
    _rendertime!(imagebuffer, textconfig(o), currenttime(data))
    return imagebuffer
end

_valn(::Nothing, n) = nothing
_valn(vals, n) = vals[n]

@noinline _grididnotinkeyserror(grid_id, grids) =
    throw(ArgumentError("$grid_id is not in $(keys(grids))"))
@noinline _wronglengtherror(f, ngrids, n) =
    throw(ArgumentError("Number of grids ($ngrids) and length of $f ($n) must be the same"))


# Automatically choose an image renderer

function autorenderer(init, scheme)
    Image(first(_iterableschemes(scheme)), ZEROCOL, MASKCOL)
end
function autorenderer(init::NamedTuple, scheme)
    rows = length(init) ÷ 4 + 1
    cols = (length(init) - 1) ÷ rows + 1
    layout = reshape([keys(init)...], (rows, cols))
    renderer = autorenderer.(values(init), _iterableschemes(scheme))
    Layout(layout, renderer)
end

_iterableschemes(::Nothing) = (ObjectScheme(),)
_iterableschemes(schemes::Union{Tuple,NamedTuple,AbstractArray}) = schemes
_iterableschemes(scheme) = (scheme,)

# Coll conversion tools

# Set a value to be between zero and one, before converting to Color.
# min and max of `nothing` are assumed to be 0 and 1.
normalise(x::X, minv, maxv) where X = max(min((x - minv) / (maxv - minv), oneunit(X)), zero(X))
normalise(x::X, minv, maxv::Nothing) where X = max((x - minv) / (oneunit(X) - minv), zero(X))
normalise(x::X, minv::Nothing, maxv) where X = min(x / maxv, oneunit(X), oneunit(X))
normalise(x, minv::Nothing, maxv::Nothing) = x
normalise(x; min=Nothing, max=Nothing) = normalise(x, min, max)

# Rescale a value between 0 and 1 to be between `min` and `max`.
# This can be used to shrink the range of a colorsheme that is displayed.
# min and max of `nothing` are assumed to be 0 and 1.
scale(x, min, max) = x * (max - min) + min
scale(x, ::Nothing, max) = x * max
scale(x, min, ::Nothing) = x * (oneunit(typeof(min)) - min) + min
scale(x, ::Nothing, ::Nothing) = x

to_rgb(vals::Tuple) = ARGB32(vals...)
to_rgb(val::Real) = ARGB32(RGB(val))
to_rgb(val::Color) = ARGB32(val)
to_rgb(val::ARGB32) = val
# Handle external colorshemes, such as from ColorSchemes.jl
to_rgb(scheme, val::Real) = to_rgb(get(scheme, val))