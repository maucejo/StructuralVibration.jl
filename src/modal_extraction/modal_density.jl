"""
    FDFit(window)

Finite difference fit for modal density estimation

**Field**
* `window::Int`: Parameter is the number of modes to consider in the finite difference approximation (default is 5)
"""
@show_data struct FDFit
    window ::Int

    FDFit(window = 5) = new(window)
end

"""
    PolyFit()

Third-order polynomial fit for modal density estimation
"""
struct PolyFit end

"""
    Lowess(window)

Lowess (Local Weighted Scatterplot Smoothing) method for modal density estimation

**Field**
* `window::Float64`: Fraction of modes to consider in the local regression (default is 0.3)
"""
@show_data struct Lowess
    window::Float64

    Lowess(window = 0.3) = new(window)
end

"""
    CSFit()

Cubic spline interpolation method for modal density estimation
"""
@show_data struct CSFit end

"""
    modal_density(frequencies, method::FDFit)
    modal_density(frequencies, method::PolyFit)
    modal_density(frequencies, method::Lowess)
    modal_density(frequencies, method::CSFit)

Compute the modal density of a structure from the knowledge of its resonance frequencies using different methods

**Inputs**
* `frequencies::Vector{Real}`: Vector of resonance frequencies
* `method`: Method to use for modal density estimation
    * `FDFit`: Finite difference fitting method
    * `PolyFit`: Third order polynomial global fitting method
    * `Lowess`: Lowess method
    * `CSFit`: Cubic spline interpolation method

**Outputs**
* `md::Vector{Real}`: Modal density values corresponding to the input frequencies
* `freqs::Vector{Real}`: Sorted frequencies corresponding to the modal density values
"""
function modal_density(frequencies::Vector{T}, method::FDFit) where T <: Real
    # Validation
    length(frequencies) < 2 && throw(ArgumentError("At least 2 frequencies required"))
    any(frequencies .<= 0) && throw(ArgumentError("All frequencies must be positive"))

    # Sort frequencies
    window_size = method.window
    freqs = sort(frequencies)
    n_modes = length(freqs)

    md = zero.(freqs)
    for i in eachindex(freqs)
        # Local window bounds
        i_start = max(1, i - window_size)
        i_end = min(n_modes, i + window_size)

        # Compute local derivative
        ΔN = i_end - i_start
        Δf = freqs[i_end] - freqs[i_start]

        md[i] = ΔN / Δf
    end

    return md, freqs
end

function modal_density(frequencies::Vector{T}, method::PolyFit) where T <: Real
    # Validation
    length(frequencies) < 2 && throw(ArgumentError("At least 2 frequencies required"))
    any(frequencies .<= 0) && throw(ArgumentError("All frequencies must be positive"))

    # Sort frequencies
    freqs = sort(frequencies)
    n_modes = length(freqs)
    modes_idx = 1:n_modes

    # Polynomial fitting approach (3rd order)
    t = (freqs .- freqs[1]) / (freqs[end] - freqs[1])  # Normalize to [0, 1]

    # Fit polynomial: N(t) = a₀ + a₁t + a₂t² + a₃t³
    A = [one.(freqs) t t.^2 t.^3]
    coeffs = A \ modes_idx

    # Derivative: dN/dt = a₁ + 2a₂t + 3a₃t²
    # Derivative: dt/df = 1/(f_max - f_min)
    # dN/df = dN/dt * dt/df
    dN_dt = coeffs[2] .+ 2coeffs[3]*t .+ 3coeffs[4]*t.^2
    df_dt = freqs[end] - freqs[1]
    md = dN_dt ./ df_dt

    return md, freqs
end

function modal_density(frequencies::Vector{T}, method::Lowess) where T <: Real
    # Validation
    length(frequencies) < 2 && throw(ArgumentError("At least 2 frequencies required"))
    any(frequencies .<= 0) && throw(ArgumentError("All frequencies must be positive"))

    # Sort frequencies
    freqs = sort(frequencies)
    n_modes = length(freqs)
    modes_idx = 1:n_modes

    window_frac = method.window
    window_size = max(2, Int(ceil(n_modes * window_frac)))

    md = zero.(freqs)
    for (i, fi) in enumerate(freqs)
        # Local window bounds
        i_start = max(1, i - window_size÷2)
        i_end = min(n_modes, i_start + window_size - 1)

        # Adjust start index to center the window around i
        i_start = max(1, i_end - window_size + 1)

        # Tricube weighting function
        idx = i_start:i_end
        t = abs.(idx .- i) / (i_end - i_start + 1) # Normalize t to [0, 1]
        weights = (1 .- t.^3).^3
        weights[t .> 1] .= 0.

        # Weighted linear regression to estimate local slope (modal density)
        A = [ones(length(idx)) (freqs[idx] .- fi)]
        W = Diagonal(weights)
        coeffs = (A'W * A) \ (A'W * modes_idx[idx])

        if freqs[i_end] ≠ freqs[i_start]
            md[i] = coeffs[2]
        end
    end

    return md, freqs
end

function modal_density(frequencies::Vector{T}, method::CSFit) where T <: Real
    # Validation
    length(frequencies) < 2 && throw(ArgumentError("At least 2 frequencies required"))
    any(frequencies .<= 0) && throw(ArgumentError("All frequencies must be positive"))

    # Sort frequencies
    freqs = sort(frequencies)
    n_modes = length(freqs)
    modes_idx = 1:n_modes

    # Interpolation approach
    itp = CubicSpline(modes_idx, freqs)
    md = DataInterpolations.derivative.(Ref(itp), freqs) # Compute derivative (modal density)

    return md, freqs
end