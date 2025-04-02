"""
    gradient(f, t)

Compute the gradient of a function `f` at points `t`.

**Inputs**
- `f::Vector{Real}`: Function values
- `t`: Points at which to evaluate the gradient
- `method`: Interpolation method
    - `:linear`: Linear interpolation
    - `:cubic`: Cubic spline interpolation (default)

**Output**
- `df`: Gradient of the vector `f` at points `t`

**Note**
- If `method` is `:cubic`, `t` must be an `AbstractRange`.
"""
function gradient(f::Vector{T}, t; method = :cubic) where {T <: Real}

    if (method == :cubic) && (t isa AbstractRange)
        itp = cubic_spline_interpolation(t, f)
    elseif method == :linear
        itp = linear_interpolation(t, f)
    end

    return only.(Interpolations.gradient.(Ref(itp), t))
end

"""
    gradient(f, t) where {T <: Real}

Compute the gradient of a function `f` at points `t`.

**Inputs**
- `f::Matrix{Real}`: Function values
- `t`: Points at which to evaluate the gradient
- `method`: Interpolation method
    - `:linear`: Linear interpolation
    - `:cubic`: Cubic spline interpolation (default)

**Output**
- `df`: Gradient of the each row of the matrix `f` at points `t`

**Note**
- If `method` is `:cubic`, `t` must be an `AbstractRange`.
"""
function gradient(f::Matrix{T}, t; method = :cubic) where {T <: Real}
    nx, nt = size(f)
    df = zeros(nx, nt)
    for i in 1:nx
        df[i, :] = gradient(f[i, :], t, method = method)
    end

    return df
end

"""
    curvature(x, y, p) where {T <: Real}

Compute the curvature of a parametric curve `(x(t), y(t))` at points `t`.

**Inputs**
- `x`: x-coordinates of the curve
- `y`: y-coordinates of the curve
- `p`: Points at which to evaluate the curvature
- `method`: Interpolation method
    - `:linear`: Linear interpolation
    - `:cubic`: Cubic spline interpolation (default)

**Output**
- `curvature`: Curvature of the curve at points `p`

**Note**
- If `method` is `:cubic`, `t` must be an `AbstractRange`.
"""
function curvature(x::Vector{T}, y::Vector{T}, p; method = :cubic) where {T <: Real}
    dx = gradient(x, p, method = method)
    ddx = gradient(dx, p, method = method)

    dy = gradient(y, p, method = method)
    ddy = gradient(dy, p, method = method)

    return (dx.*ddy - dy.*ddx)./(dx.^2 + dy.^2).^(3/2)
end