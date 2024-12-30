function gradient(f; h = 1.)
    n = length(f)
    df = zeros(n)

    # First point
    df[1] = (f[2] - f[1])/h

    # Last point
    df[end] = (f[end] - f[end-1])/h

    # Central difference
    df[2:end-1] = (f[3:end] - f[1:end-2])/2h

    # Last point
    return df
end

function gradient(f, t)
    itp = LinearInterpolation(t, f)

    return only.(Interpolations.gradient.(Ref(itp), t))
end