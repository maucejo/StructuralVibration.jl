"""
    detrend(t, y, order = 1, bp = [])

Detrend a signal `y` with respect to time `t` using a polynomial of order `order`.

**Inputs**
* `t`: Time vector
* `y`: Signal to be detrended
* `order`: Order of the polynomial (default is 1)
* `bp`: Breakpoints for the polynomial (default is empty)

**Output**
- `y_detrended`: Detrended signal

**Notes**

`order` can be a vector if the trend is different over the segments defined by the breakpoints `bp`.

The breakpoints `bp` are the points where the polynomial order changes. If `bp` is empty, a single polynomial is fitted to the entire signal.
"""
function detrend(t, y, order = 1, bp = eltype(y)[])
    if length(bp) == 0
        return y .- polyval(polyfit(t, y, order[1]), t)
    else
        idb = findall(@. bp ≥ t[1] || bp ≤ t[end])

        # If there is only one breakpoint, convert it to a scalar
        length(idb) == 1 ? idb = only(idb) : nothing

        bpi = unique([t[1]; bp[idb]; t[end]])

        if length(order) == 1
            order = order*ones(length(bpi) - 1)
        else
            if length(order) != length(bpi) - 1
                throw(ArgumentError("The number of orders must be equal to the number of segments."))
            end
        end

        y_detrended = similar(y)
        for (i, ordi) in enumerate(order)
            idx = findall(@. bpi[i] ≤ t ≤ bpi[i + 1])

            y_detrended[idx] .= y[idx] .- polyval(polyfit(t[idx], y[idx], ordi), t[idx])
        end

        return y_detrended
    end
end

"""
    polyfit(x, y, order = 1)

Fit a polynomial of order `order` to the data `x` and `y`.

**Inputs**
* `x`: Independent variable
* `y`: Dependent variable
* `order::Real`: Order of the polynomial (default is 1)

**Output**
* `p`: Coefficients of the polynomial
"""
function polyfit(x, y, order::Int = 1)
    if order == 0
        return [mean(y)]
    else
        orders = 0:order
        A = [xi^ordj for xi in x, ordj in reverse(orders)]

        return (A'A)\(A'y)
    end
end

"""
    polyval(p, x)

    Evaluate a polynomial of coefficients `p` at `x`

**Inputs**
* `p::Vector{Real}`: Coefficients of the polynomial
* `x`: Value at which the polynomial is evaluated

**Output**
* `y`: Value of the polynomial at `x`
"""
function polyval(p, x)
    y = zeros(eltype(p), length(x))
    orders = 0:length(p) - 1

    for (ci, ordi) in zip(p, reverse(orders))
        @. y += ci*x^ordi
    end

    return y
end