"""
    gradient(f::Vector{Real}, t; method = :cubic)
    gradient(f::Matrix{Real}, t; method = :cubic, dims = 1)

Compute the gradient of a function `f` at points `t`.

**Inputs**
* `f`: Function values
* `t`: Points at which to evaluate the gradient
* `method`: Interpolation method
    - `:linear`: Linear interpolation
    - `:cubic`: Cubic spline interpolation (default)
* `dims`: Dimension along which to compute the gradient
    * `1`: Rows (default)
    * `2`: Columns

**Output**
* `df`: Gradient of the vector `f` at points `t`

**Note**

If `method` is `:cubic`, `t` must be an `AbstractRange`.
"""
function gradient(f::Vector{T}, t; method = :cubic) where {T <: Real}
    if (method == :cubic)
        itp = CubicSpline(f, t)
    elseif method == :linear
        itp = LinearInterpolation(f, t)
    end

    return DataInterpolations.derivative.(Ref(itp), t)
end

function gradient(f::Matrix{T}, t; method = :cubic, dims = 1) where {T <: Real}
    nx, nt = size(f)
    df = zeros(nx, nt)

    if dims == 1
        for i in 1:nx
            df[i, :] = gradient(f[i, :], t, method = method)
        end
    else
        for i in 1:nt
            df[:, i] = gradient(f[:, i], t, method = method)
        end
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