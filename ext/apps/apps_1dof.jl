"""
    free_response_1dof(m = 1., t = 0.:1e-3:2.)

Launch an interactive application to visualize the free response of a single-degree-of-freedom (SDOF) system

**Inputs**
- `m`: Mass of the SDOF system (default: 1.0)
- `t`: Time vector over which the response is computed (default: 0.0
    to 2.0 seconds with a step of 1e-3)

**Output**
- App
"""
function free_response_1dof(m = 1., t = 0.:1e-3:2.)
    function response(f0, ξ, x0, v0)
        prob = SdofFreeTimeProblem(Sdof(m, f0, ξ), [x0, v0], t)
        return solve(prob).u
    end

    ## Application
    GLMakie.activate!(title = "Free response of a SDOF system")
    theme = theme_choice(:makie, fonts = (regular = "Lato", bold = "Lato Bold"))

    fig = with_theme(theme) do
        fig = Figure(size = (1280, 720), figure_padding = 30.)

        Label(fig[1, 1:2], "Free response of a SDOF system", fontsize = 38., halign = :left, tellwidth = false, color = :gray20)

        # Graph
        ax = Axis(fig[2, 1], xlabel = "Time (s)", ylabel = "Displacement (m)", yticklabelspace = 60.,
            xlabelsize = 20., ylabelsize = 20., xticklabelsize = 15., yticklabelsize = 15.,
            topspinevisible = false, rightspinevisible = false)

        # Control panel
        frame = GridLayout(fig[2, 2], tellheight = false)
        Box(frame[1, 1], color = :transparent, strokewidth = 1., strokecolor = :gray30, cornerradius = 6.)
        panel = GridLayout(frame[1, 1], alignmode = Outside(18.))

        f0 = add_param!(panel, 1, 1.:1.:50., 10.,
            v -> "Natural frequency : $(round(Int, v)) Hz",
            [0., 9., 19., 29., 39., 49.]/49., ["1", "10", "20", "30", "40", "50"])

        ξ = add_param!(panel, 2, 0.:0.01:1.5, 0.,
            v -> "Damping ratio : $v",
            [0., 0.5, 1., 1.5]/1.5, ["0", "0.5", "1", "1.5"])

        x0 = add_param!(panel, 3, -1.:0.1:1., 0.5,
            v -> "Initial displacement : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"])

        v0 = add_param!(panel, 4, -100.:1.:100., 0.,
            v -> "Initial velocity : $(round(Int, v)) m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-100", "-50", "0", "50", "100"])

        # Checkbox to freeze the y-axis scale
        option = GridLayout(panel[9, 1], halign = :left)
        checkbox = Checkbox(option[1, 1], checked = false, size = 18., checkboxcolor_checked = colorant"#7fc4e8",
            checkboxstrokecolor_checked = colorant"#7fc4e8", checkboxstrokecolor_unchecked = :gray40)
        Label(option[1, 2], "Freeze y-axis", fontsize = 16., color = :gray20)
        colgap!(option, 8.)
        y_fixed = checkbox.checked

        for i in 1:8
            rowgap!(panel, i, isodd(i) ? 8. : 38.)
        end

        colsize!(fig.layout, 2, Relative(0.26))

        # Response
        x = lift(response, f0, ξ, x0, v0)

        hlines!(ax, 0., color = :gray30, linewidth = 1.)
        lines!(ax, t, x, linewidth = 2.)
        xlims!(ax, t[begin], t[end])

        # Automatic scaling, unless the scale is fixed (rescaled on the curve when released)
        onany(x, y_fixed, update = true) do xs, fixed
            fixed && return
            xmin, xmax = extrema(xs)
            δ = max(0.05*(xmax - xmin), 1e-3)
            ylims!(ax, xmin - δ, xmax + δ)
        end

        fig
    end

    screen = display(fig)

    # Keep the window open when the script is launched from a terminal
    if abspath(PROGRAM_FILE) == @__FILE__
        wait(screen)
    end
end


"""
    forced_response_1dof(m, t, r, Fmax, xbmax)

Launch an interactive application to visualize the forced response of a single-degree-of-freedom (SDOF) system

**Inputs**
- `m`: Mass of the system (default: 1.0)
- `t`: Time vector (default: 0.:1e-3:2.)
- `r`: Reduced frequency vector (default: 0.:0.005:5.)
- `Fmax`: Maximum excitation force (default: 20.0)
- `xbmax`: Maximum base displacement (default: 2.0)

**Output**
- App
"""
function forced_response_1dof(m = 1., t = 0.:1e-3:2., r = 0.:0.005:5., Fmax = 20., xbmax = 2.)

    # Excitation amplitude (a : relative position on the slider)
    amplitude(a, type_exc) = type_exc == :force ? round(Fmax*a) : round(xbmax*a, digits = 1)

    # Forced response of a SDOF system
    function response(f0, ξ, f, A, x0, v0, type_exc)
        prob = SdofHarmonicTimeProblem(Sdof(m, f0, ξ), [x0, v0], t, A, f, type_exc)
        return solve(prob).u
    end

    # Transfer/Transmissibility function of a SDOF system
    function transfer(f0, ξ, type_exc)
        prob = SdofFRFProblem(Sdof(m, f0, ξ), f0*r, type_exc = type_exc)
        return solve(prob).u
    end

    GLMakie.activate!(title = "Forced response of a SDOF system")
    theme = theme_choice(:makie, fonts = (regular = "Lato", bold = "Lato Bold"))

    fig = with_theme(theme) do
        fig = Figure(size = (1400, 840), figure_padding = 30.)

        Label(fig[1, 1:2], "Forced response of a SDOF system", fontsize = 34., halign = :left, tellwidth = false, color = :gray20)

        axis_style = (xlabelsize = 18., ylabelsize = 18., xticklabelsize = 14., yticklabelsize = 14.,
            titlesize = 18., titlealign = :left, titlefont = :regular,
            topspinevisible = false, rightspinevisible = false)

        # Graphs
        ax_x = Axis(fig[2, 1]; title = "Forced response", xlabel = "Time (s)", ylabel = "Displacement (m)", yticklabelspace = 55., axis_style...)
        ax_amp = Axis(fig[3, 1]; xlabel = "Reduced frequency f/f₀", yscale = log10, yticklabelspace = 55.,
            ytickformat = vs -> [rich("10", superscript(string(round(Int, log10(v))))) for v in vs],
            yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9), axis_style...)
        ax_phi = Axis(fig[3, 2]; xlabel = "Reduced frequency f/f₀", ylabel = "Phase (rad)",
            yticks = ([-π, -π/2, 0., π/2], ["-π", "-π/2", "0", "π/2"]), axis_style...)

        # Control panel
        frame = GridLayout(fig[2, 2], tellheight = false)
        Box(frame[1, 1], color = :transparent, strokewidth = 1., strokecolor = :gray30, cornerradius = 6.)
        panel = GridLayout(frame[1, 1], alignmode = Outside(18.))

        Label(panel[1, 1], "Excitation type", halign = :left, fontsize = 16., color = :gray20)
        menu = Menu(panel[1, 2], options = [("External force", :force), ("Base motion", :base)], default = "External force", fontsize = 15.)
        type_exc = menu.selection

        a = add_param!(panel, 1, 1, 0.:0.05:1., 0.5,
            (a, e) -> e == :force ? "Force amplitude : $(round(Int, amplitude(a, e))) N" : "Base amplitude : $(amplitude(a, e)) m",
            [0., 0.25, 0.5, 0.75, 1.], lift(e -> e == :force ? ["0", "5", "10", "15", "20"] : ["0", "0.5", "1", "1.5", "2"], type_exc),
            deps = (type_exc,), offset = 1)

        f = add_param!(panel, 1, 2, 1.:1.:50., 5.,
            v -> "Excitation frequency : $(round(Int, v)) Hz",
            [0., 9., 19., 29., 39., 49.]/49., ["1", "10", "20", "30", "40", "50"], offset = 1)

        f0 = add_param!(panel, 2, 1, 1.:1.:50., 10.,
            v -> "Natural frequency : $(round(Int, v)) Hz",
            [0., 9., 19., 29., 39., 49.]/49., ["1", "10", "20", "30", "40", "50"], offset = 1)

        ξ = add_param!(panel, 2, 2, 0.:0.01:1.5, 0.1,
            v -> "Damping ratio : $v",
            [0., 0.5, 1., 1.5]/1.5, ["0", "0.5", "1", "1.5"], offset = 1)

        x0 = add_param!(panel, 3, 1, -1.:0.1:1., 0.5,
            v -> "Initial displacement : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"], offset = 1)

        v0 = add_param!(panel, 3, 2, -100.:1.:100., 0.,
            v -> "Initial velocity : $(round(Int, v)) m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-100", "-50", "0", "50", "100"], offset = 1)

        # Checkbox to freeze the y-axis scale
        option = GridLayout(panel[8, 1:2], halign = :left)
        checkbox = Checkbox(option[1, 1], checked = false, size = 18., checkboxcolor_checked = colorant"#7fc4e8",
            checkboxstrokecolor_checked = colorant"#7fc4e8", checkboxstrokecolor_unchecked = :gray40)
        Label(option[1, 2], "Freeze y-axis", fontsize = 16., color = :gray20)
        colgap!(option, 8.)
        y_fixed = checkbox.checked

        for i in 1:6
            rowgap!(panel, i, isodd(i) ? 14. : 8.)
        end
        for i in 3:2:7
            rowgap!(panel, i, 34.)
        end
        colgap!(panel, 1, 30.)

        colsize!(fig.layout, 1, Relative(0.5))
        colsize!(fig.layout, 2, Relative(0.5))

        # Forced response
        A = lift(amplitude, a, type_exc)
        x = lift(response, f0, ξ, f, A, x0, v0, type_exc)

        hlines!(ax_x, 0., color = :gray30, linewidth = 1.)
        lines!(ax_x, t, x, linewidth = 2.)
        xlims!(ax_x, t[begin], t[end])

        # Automatic scaling, unless the scale is fixed (rescaled on the curve when released)
        onany(x, y_fixed, update = true) do xs, fixed
            fixed && return
            xmin, xmax = extrema(xs)
            δ = max(0.05*(xmax - xmin), 1e-3)
            ylims!(ax_x, xmin - δ, xmax + δ)
        end

        # Transfer/Transmissibility function
        H = lift(transfer, f0, ξ, type_exc)
        amp = lift(Hs -> [isfinite(z) ? abs(z) : NaN for z in Hs], H)
        phi = lift(Hs -> [isfinite(z) ? angle(z) : NaN for z in Hs], H)

        # Locate the excitation frequency on the graphs
        r_exc = lift((f, f0) -> [f/f0], f, f0)

        connect!(ax_amp.title, lift(e -> e == :force ? "Transfer function H" : "Transmissibility function T", type_exc))
        connect!(ax_phi.title, lift(e -> e == :force ? "Phase of H" : "Phase of T", type_exc))
        connect!(ax_amp.ylabel, lift(e -> e == :force ? "Amplitude (m/N)" : "Amplitude", type_exc))

        for ax in (ax_amp, ax_phi)
            vlines!(ax, r_exc, color = colorant"#e36837", linewidth = 1.5, linestyle = :dash)
            xlims!(ax, r[begin], r[end])
        end

        lines!(ax_amp, r, amp, linewidth = 2.)
        lines!(ax_phi, r, phi, linewidth = 2.)
        ylims!(ax_phi, -π - 0.15, π/2 + 0.15)

        onany(amp, y_fixed, update = true) do as, fixed
            fixed && return
            amin, amax = extrema(filter(isfinite, as))
            kmin, kmax = floor(Int, log10(amin)), ceil(Int, log10(amax))
            ax_amp.yticks = exp10.(kmin:kmax)
            ylims!(ax_amp, exp10(kmin), exp10(kmax))
        end

        fig
    end

    screen = display(fig)

    # Keep the window open when the script is run from a terminal
    if abspath(PROGRAM_FILE) == @__FILE__
        wait(screen)
    end
end