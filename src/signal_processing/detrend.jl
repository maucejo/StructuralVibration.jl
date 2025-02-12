"""
detrend(t, y, order=1)

Detrend a signal `y` with respect to time `t` using a polynomial of order `order`.
"""
function detrend(t::AbstractArray, y::AbstractArray, order = 1, bp = Float64[])
    if length(bp) == 0
        return y - polyval(polyfit(t, y, order[1]), t)
    else
        idb = findall(@. bp ≥ t[1] || bp ≤ t[end])

        # If there is only one breakpoint, convert it to a scalar
        length(idb) == 1 ? idb = only(idb) : nothing

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
    polyfit(x, y, order = 1)

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
    else
        orders = 0:order
        A = [xi^ordj for xi in x, ordj in reverse(orders)]

        return (A'A)\(A'y)
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
    y = zeros(length(x))
    orders = 0:length(p) - 1

    for (ci, ordi) in zip(p, reverse(orders))
        @. y += ci*x^ordi
    end

    return y
end