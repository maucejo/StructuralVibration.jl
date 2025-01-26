"""
    gradient(f::Vector{Float64}, t)

Compute the gradient of a function `f` at points `t`.

# Inputs
- `f::Vector{Float64}`: Function values
- `t`: Points at which to evaluate the gradient

# Output
- `df`: Gradient of the vector `f` at points `t`
"""
function gradient(f::Vector{Float64}, t)

    itp = linear_interpolation(t, f)

    return only.(Interpolations.gradient.(Ref(itp), t))
end

"""
    gradient(f::Matrix{Float64}, t)

Compute the gradient of a function `f` at points `t`.

# Inputs
- `f::Matrix{Float64}`: Function values
- `t`: Points at which to evaluate the gradient

# Output

- `df`: Gradient of the each row of the matrix `f` at points `t`
"""
function gradient(f::Matrix{Float64}, t)
    nx, nt = size(f)
    df = zeros(nx, nt)
    for i in 1:nx
        df[i, :] = gradient(f[i, :], t)
    end

    return df
end

function curvature(x::Vector{Float64}, y::Vector{Float64}, p)
    dx = gradient(x, p)
    ddx = gradient(dx, p)

    dy = gradient(y, p)
    ddy = gradient(dy, p)

    return (dx.*ddy - dy.*ddx)./(dx.^2 + dy.^2).^(3/2)
end