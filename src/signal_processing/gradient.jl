"""
    gradient(f, t)

Compute the gradient of a function `f` at points `t`.

# Inputs
- `f::Vector{Real}`: Function values
- `t::AbstractRange`: Points at which to evaluate the gradient

# Output
- `df`: Gradient of the vector `f` at points `t`
"""
function gradient(f::Vector{T}, t::AbstractRange) where {T <: Real}

    itp = cubic_spline_interpolation(t, f)

    return only.(Interpolations.gradient.(Ref(itp), t))
end

"""
    gradient(f, t) where {T <: Real}

Compute the gradient of a function `f` at points `t`.

# Inputs
- `f::Matrix{Real}`: Function values
- `t::AbstractRange`: Points at which to evaluate the gradient

# Output

- `df`: Gradient of the each row of the matrix `f` at points `t`
"""
function gradient(f::Matrix{T}, t::AbstractRange) where {T <: Real}
    nx, nt = size(f)
    df = zeros(nx, nt)
    for i in 1:nx
        df[i, :] = gradient(f[i, :], t)
    end

    return df
end

"""
    curvature(x, y, p) where {T <: Real}

Compute the curvature of a parametric curve `(x(t), y(t))` at points `t`.

# Inputs
- `x`: x-coordinates of the curve
- `y`: y-coordinates of the curve
- `p`: Points at which to evaluate the curvature

# Output
- `curvature`: Curvature of the curve at points `p`
"""
function curvature(x::Vector{T}, y::Vector{T}, p) where {T <: Real}
    dx = gradient(x, p)
    ddx = gradient(dx, p)

    dy = gradient(y, p)
    ddy = gradient(dy, p)

    return (dx.*ddy - dy.*ddx)./(dx.^2 + dy.^2).^(3/2)
end