# See interface docs
@inline inbounds(data::Union{GridData,AbstractSimData}, I::Tuple) = inbounds(data, I...)
@inline inbounds(data::Union{GridData,AbstractSimData}, I...) = 
    _inbounds(boundary(data), gridsize(data), I...)

@inline function _inbounds(boundary::BoundaryCondition, size::Tuple, i1, i2)
    a, inbounds_a = _inbounds(boundary, size[1], i1)
    b, inbounds_b = _inbounds(boundary, size[2], i2)
    (a, b), inbounds_a & inbounds_b
end
@inline _inbounds(::Remove, size::Number, i::Number) = i, _isinbounds(size, i)
@inline function _inbounds(::Wrap, size::Number, i::Number)
    if i < oneunit(i)
        size + rem(i, size), true
    elseif i > size
        rem(i, size), true
    else
        i, true
    end
end
@inline function DG._inbounds(::Reflect, size::Number, i::Number)
    if i < oneunit(i)
        return 2 - i, true
    elseif i > size
        return 2*size - i, true
    else
        return i, true
    end
end

# See interface docs
@inline isinbounds(data::Union{GridData,AbstractSimData}, I::Tuple) = isinbounds(data, I...)
@inline isinbounds(data::Union{GridData,AbstractSimData}, I...) = _isinbounds(gridsize(data), I...)

@inline _isinbounds(size::Tuple, I...) = all(map(_isinbounds, size, I))
@inline _isinbounds(size, i) = i >= one(i) && i <= size