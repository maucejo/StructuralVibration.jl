using StructuralVibration, Statistics
@usingany CairoMakie, PrettyTables

## Clamped-clamped string
# Geometry and material parameters
L = 5.        # Length
b = 0.03      # Width
h = 0.01      # Thickness
S = b*h       # Cross-section area
ρ = 2698.   # Density

# Create a Strings struct
strings = Strings(L, S, 1e4, ρ)

# Reference modal density
md_ref_string = round(modal_density(strings), digits = 5)

# Modal density estimation from the frequencies
ωn_string = modefreq(strings, 1e3)[1]
fn_string = ωn_string./2π

md_fd_string = round.(modal_density(fn_string, FDFit())[1], digits = 5)
md_poly_string = round.(modal_density(fn_string, PolyFit())[1], digits = 5)
md_lowess_string = round.(modal_density(fn_string, Lowess())[1], digits = 5)
md_interp_string = round.(modal_density(fn_string, CSFit())[1], digits = 5)

## Display results
fig_string = begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Frequency (Hz)", ylabel = "Modal density (modes/Hz)", title = "Clamped-Clamped string")
    sc1 = scatter!(ax, fn_string, md_fd_string)
    sc2 = scatter!(ax, fn_string, md_poly_string, marker = :star5)
    sc3 = scatter!(ax, fn_string, md_lowess_string, marker = :xcross)
    sc4 = scatter!(ax, fn_string, md_interp_string, marker = :diamond)
    h1 = hlines!(ax, [md_ref_string], color = :black)
    xlims!(ax, minimum(fn_string), maximum(fn_string))
    axislegend(ax, [h1, sc1, sc2, sc3, sc4], ["Reference", "FDFit", "PolyFit", "Lowess", "Interp"])

    fig
end

# Table Construction
table_string = begin
    data = [md_ref_string mean(md_fd_string) mean(md_poly_string) mean(md_lowess_string) mean(md_interp_string)]
    header = ["Reference", "FDFit", "PolyFit", "Lowess", "CSFit"]
    pretty_table(data, column_labels = header, alignment = :c)
end

## Simply-supported beam
# Geometry and material parameters
L = 1.        # Length
b = 0.03      # Width
h = 0.01      # Thickness
S = b*h       # Cross-section area
Iz = b*h^3/12 # Moment of inertia
E = 2.1e11  # Young's modulus
ρ = 7850.   # Density

# Create a Beam struct
beam = Beam(L, S, Iz, E, ρ)
fmax = 5000.

# Theoretical resonance frequencies
ωn_beam = modefreq(beam, fmax)[1]
fn_beam = ωn_beam./2π

# Reference modal density
freq_beam = LinRange(minimum(fn_beam), maximum(fn_beam), 1000)
md_ref_beam = round.(modal_density.(Ref(beam), freq_beam), digits = 4)

# Modal density estimation from the frequencies
ωn_beam = modefreq(beam, fmax)[1]
fn_beam = ωn_beam./2π

md_fd_beam = round.(modal_density(fn_beam, FDFit(2))[1], digits = 4)
md_poly_beam = round.(modal_density(fn_beam, PolyFit())[1], digits = 4)
md_lowess_beam = round.(modal_density(fn_beam, Lowess(0.1))[1], digits = 4)
md_interp_beam = round.(modal_density(fn_beam, CSFit())[1], digits = 4)

fig_beam = begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Frequency (Hz)", ylabel = "Modal density (modes/Hz)", title = "Simply-supported beam")
    lines!(ax, freq_beam, md_ref_beam, color = :black, label = "Reference")
    scatter!(ax, fn_beam, md_fd_beam, label = "FDFit")
    scatter!(ax, fn_beam, md_poly_beam, marker = :star5, label = "PolyFit")
    scatter!(ax, fn_beam, md_lowess_beam, marker = :xcross, label = "Lowess")
    scatter!(ax, fn_beam, md_interp_beam, marker = :diamond, label = "CSFit")
    xlims!(ax, minimum(fn_beam), maximum(fn_beam))
    axislegend(ax)

    fig
end

md_ref_b = round.(modal_density.(Ref(beam), fn_beam), digits = 4)
table_beam = begin
    data = [round(mean(md_ref_b), digits = 4) round(mean(md_fd_beam), digits = 4) round(mean(md_poly_beam), digits = 4) round(mean(md_lowess_beam), digits = 4) round(mean(md_interp_beam), digits = 4)]
    header = ["Reference", "FDFit", "PolyFit", "Lowess", "CSFit"]
    pretty_table(data, column_labels = header, alignment = :c)
end