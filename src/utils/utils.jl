# Convenient function to create arrays of undefs
"""
    undefs(::Type{T}, dims::Int...) where T
    undefs(dims::Int...)

Create an array of undefs of type `T` and dimensions `dims`. If `T` is not provided, the default type is `Float64`.

# Inputs
- `T::Type`: Type of the array
- `dims::Int...`: Dimensions of the array

# Output
- `arr`: Array of undefs
"""
undefs(::Type{T}, dims::Int...) where T = Array{T}(undef, dims)
undefs(dims::Int...) = Array{Float64}(undef, dims)

"""
detrend(t, y, order=1)

Detrend a signal `y` with respect to time `t` using a polynomial of order `order`.
"""
function detrend(t::AbstractArray, y::AbstractArray, order = 1, bp = Float64[])
    if length(bp) == 0
        return y - polyval(polyfit(t, y, order[1]), t)
    else
        idb = findall(@. bp ≥ t[1] || bp ≤ t[end])
        if length(idb) == 1
            idb = only(idb)
        end
        bpi = unique([t[1]; bp[idb]; t[end]])

        if length(order) == 1
            order = order*ones(length(bpi) - 1)
        end

        y_detrended = undefs(length(y))

        for (i, ordi) in enumerate(order)
            idx = findall(@. bpi[i] ≤ t ≤ bpi[i + 1])

            y_detrended[idx] = y[idx] - polyval(polyfit(t[idx], y[idx], ordi), t[idx])
        end

        return y_detrended
    end
end

"""
    polyfit(x, y, order=1)

Fit a polynomial of order `order` to the data `x` and `y`.

# Inputs
- `x`: Independent variable
- `y`: Dependent variable
- `order`: Order of the polynomial (default is 1)

# Output
- `p`: Coefficients of the polynomial
"""
function polyfit(x::AbstractArray, y::AbstractArray, order::Int = 1)
    if order == 0
        return mean(y)

    elseif order == 1
        A = [x ones(length(x))]
        return (A'A)\(A'y)
    elseif order == 2
        A = [x.^2 x ones(length(x))]
        return (A'A)\(A'y)
    else
        error("Order > 2 not supported")
    end
end

"""
    polyval(p, x)

    Evaluate a polynomial of coefficients `p` at `x`

# Inputs
- `p`: Coefficients of the polynomial
- `x`: Value at which the polynomial is evaluated

# Output
- `y`: Value of the polynomial at `x`
"""
function polyval(p, x)
    if length(p) == 1
        return p[1]*ones(length(x))

    elseif length(p) == 2
        return @. p[1]*x + p[2]

    else
        return @. p[1]*x^2 + p[2]*x + p[3]
    end
end