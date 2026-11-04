"""
    free_response_2dof(m1, m2, t)

Launch an interactive application to visualize the free response of a two-degree-of-freedom (2DOF) system

**Inputs**
- `m1`: Mass of the first degree of freedom (default: 1.0)
- `m2`: Mass of the second degree of freedom (default: 1.0)
- `t`: Time vector over which the response is computed (default: 0.0 to 2.0 seconds with a step of 1e-3)

**Output**
- App
"""
function free_response_2dof(m1 = 1., m2 = 1., t = 0.:1e-3:2.)
    # Mass matrix
    M = mass(m1, m2)

    # Free response of a 2DOF system
    function response(f1, f2, ξ, x01, x02, v01, v02)
        prob = FreeModalTimeProblem(stiffness(m1, m2, f1, f2), M, ξ, ([x01, x02], [v01, v02]), t)
        return solve(prob).u
    end

    ## Application
    GLMakie.activate!(title="Free response of a 2DOF system")
    theme = theme_choice(:makie, fonts=(regular="Lato", bold="Lato Bold"))

    fig = with_theme(theme) do
        fig = Figure(size=(1280, 800), figure_padding=30.)

        Label(fig[1, 1:2], "Free response of a 2DOF system", fontsize=38., halign=:left, tellwidth=false, color=:gray20)

        # Graphique
        ax = Axis(fig[3, 1], xlabel="Time (s)", ylabel="Displacement (m)", yticklabelspace=60.,
            xlabelsize=20., ylabelsize=20., xticklabelsize=15., yticklabelsize=15.,
            topspinevisible=false, rightspinevisible=false)

        # Control panel
        frame = GridLayout(fig[2:3, 2], tellheight=false)
        Box(frame[1, 1], color=:transparent, strokewidth=1., strokecolor=:gray30, cornerradius=6.)
        panel = GridLayout(frame[1, 1], alignmode=Outside(18., 18., 18., 18.))

        freq_fracs = [0., 9., 19.]/19.
        freq_labels = ["1", "10", "20"]

        f1 = add_param!(panel, 1, 1.:0.5:20., 5.,
            v -> "Natural frequency of system 1 : $v Hz",
            freq_fracs, freq_labels)

        f2 = add_param!(panel, 2, 1.:0.5:20., 5.,
            v -> "Natural frequency of system 2 : $v Hz",
            freq_fracs, freq_labels)

        ξ = add_param!(panel, 3, 0.:0.005:0.5, 0.01,
            v -> "Modal damping ratio : $v",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "0.1", "0.2", "0.3", "0.4", "0.5"])

        x01 = add_param!(panel, 4, -1.:0.1:1., 0.,
            v -> "Initial displacement x₁ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"])

        x02 = add_param!(panel, 5, -1.:0.1:1., 0.5,
            v -> "Initial displacement x₂ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"])

        v01 = add_param!(panel, 6, -10.:0.5:10., 0.,
            v -> "Initial velocity v₁ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"])

        v02 = add_param!(panel, 7, -10.:0.5:10., 0.,
            v -> "Initial velocity v₂ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"])

        # Natural frequencies of the coupled system
        fn = lift((f1, f2) -> coupled_frequencies(m1, m2, f1, f2), f1, f2)
        Label(panel[15, 1], lift(fn -> "Natural frequencies : $(round(fn[1], digits = 2)) et $(round(fn[2], digits = 2)) Hz", fn),
            halign=:left, fontsize=15., color=:gray20)

        Label(panel[16, 1], "Masses : m₁ = $m1 kg et m₂ = $m2 kg", halign=:left, fontsize=15., color=:gray20)

        # Checkbox to freeze the y-axis scale
        option = GridLayout(panel[17, 1], halign=:left)
        checkbox = Checkbox(option[1, 1], checked=false, size=18., checkboxcolor_checked=colorant"#7fc4e8",
            checkboxstrokecolor_checked=colorant"#7fc4e8", checkboxstrokecolor_unchecked=:gray40)
        Label(option[1, 2], "Freeze y-axis", fontsize=16., color=:gray20)
        colgap!(option, 8.)
        y_fixed = checkbox.checked

        for i in 1:14
            rowgap!(panel, i, isodd(i) ? 8. : 34.)
        end
        rowgap!(panel, 15, 8.)
        rowgap!(panel, 16, 20.)

        colsize!(fig.layout, 2, Relative(0.3))

        # Response of the 2DOF system
        x = lift(response, f1, f2, ξ, x01, x02, v01, v02)
        x1 = lift(xs -> xs[1, :], x)
        x2 = lift(xs -> xs[2, :], x)

        hlines!(ax, 0., color=:gray30, linewidth=1.)
        lines!(ax, t, x1, color=Cycled(1), linewidth=2., label="x₁(t)")
        lines!(ax, t, x2, color=Cycled(2), linewidth=2., label="x₂(t)")
        xlims!(ax, t[begin], t[end])

        Legend(fig[2, 1], ax, orientation=:horizontal, halign=:left, framevisible=false, tellwidth=false, labelsize=16.)

        # Automatic scaling, unless the scale is fixed (rescaled on the curve when released)
        onany(x, y_fixed, update=true) do xs, fixed
            fixed && return
            xmin, xmax = extrema(xs)
            δ = max(0.05*(xmax - xmin), 1e-3)
            ylims!(ax, xmin - δ, xmax + δ)
        end

        fig
    end

    screen = display(fig)

    # Keep the window open when the script is run from a terminal
    if abspath(PROGRAM_FILE) == @__FILE__
        wait(screen)
    end
end

"""
    forced_response_2dof(m1, m2, t, freq)

Launch an interactive application to visualize the forced response of a two-degree-of-freedom (2DOF) system

**Inputs**
- `m1`: Mass of the first degree of freedom (default: 1.0)
- `m2`: Mass of the second degree of freedom (default: 1.0)
- `t`: Time vector over which the response is computed (default: 0.0 to 5.0 seconds with a step of 1e-3)
- `freq`: Frequency vector over which the transfer function is computed (default: 0.0 to 30.0 Hz with a step of 0.01)

**Output**
- App
"""
function forced_response_2dof(m1 = 1., m2 = 1., t = 0.:1e-3:5., freq = 0.:0.01:30.)
    # Mass matrix
    M = mass(m1, m2)

    # Indexes for the subscripts in the labels of the parameters
    subscripts = ["₁", "₂"]

    # Response of a 2DOF system
    function response(f1, f2, ξ, A, k, f, x01, x02, v01, v02)
        F = k == 1 ? [A, 0.] : [0., A]
        prob = HarmonicModalTimeProblem(stiffness(m1, m2, f1, f2), M, ξ, F, f, ([x01, x02], [v01, v02]), t)
        return solve(prob).u
    end

    # Transfer function of a 2DOF system
    function transfer(f1, f2, ξ)
        ωn, Φn = eigenmode(stiffness(m1, m2, f1, f2), M)
        prob = ModalFRFProblem(ωn, ξ, freq, Φn, Φn)
        return solve(prob, progress=false).u
    end

    # Phase in the range ]-2π, 0]
    function phase(z)
        isfinite(z) || return NaN
        ϕ = angle(z)
        return ϕ > 1e-9 ? ϕ - 2π : ϕ
    end

    ## Application
    GLMakie.activate!(title="Forced response of a 2DOF system")
    theme = theme_choice(:makie, fonts=(regular="Lato", bold="Lato Bold"))

    fig = with_theme(theme) do
        fig = Figure(size=(1500, 900), figure_padding=30.)

        Label(fig[1, 1:2], "Forced response of a 2DOF system", fontsize=34., halign=:left, tellwidth=false, color=:gray20)

        axis_style = (xlabelsize=18., ylabelsize=18., xticklabelsize=14., yticklabelsize=14.,
            titlesize=18., titlealign=:left, titlefont=:regular,
            topspinevisible=false, rightspinevisible=false)

        # Graphs
        time_plot = GridLayout(fig[2, 1])
        ax_x = Axis(time_plot[2, 1]; xlabel="Time (s)", ylabel="Displacement (m)", yticklabelspace=55., axis_style...)
        ax_amp = Axis(fig[3, 1]; xlabel="Frequency (Hz)", ylabel="Amplitude (dB, ref. 1 m/N)", yticklabelspace=55., axis_style...)
        ax_phi = Axis(fig[3, 2]; xlabel="Frequency (Hz)", ylabel="Phase (rad)",
            yticks=([-2π, -3π/2, -π, -π/2, 0.], ["-2π", "-3π/2", "-π", "-π/2", "0"]), axis_style...)

        # Control panel
        frame = GridLayout(fig[2, 2], tellheight=false)
        Box(frame[1, 1], color=:transparent, strokewidth=1., strokecolor=:gray30, cornerradius=6.)
        panel = GridLayout(frame[1, 1], alignmode=Outside(18.))

        Label(panel[1, 1], "Force applied on", halign = :left, tellwidth = false, fontsize = 16., color = :gray20)
        menu_k = Menu(panel[2, 1], options = [("Mass 1", 1), ("Mass 2", 2)], default = "Mass 2", fontsize = 15.)
        k = menu_k.selection


        Label(panel[1, 2], "Response observed on", halign=:left, tellwidth = false, fontsize=16., color=:gray20)
        menu_j = Menu(panel[2, 2], options=[("Mass 1", 1), ("Mass 2", 2)], default="Mass 2", fontsize=15.)
        j = menu_j.selection

        A = add_param!(panel, 1, 1, 0.:10.:500., 100.,
            v -> "Amplitude F : $(round(Int, v)) N",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "100", "200", "300", "400", "500"], offset=2)

        f = add_param!(panel, 1, 2, 0.5:0.1:30., 3.,
            v -> "Excitation freq. : $v Hz",
            [0., 9.5, 19.5, 29.5]/29.5, ["0.5", "10", "20", "30"], offset=2)

        ξ = add_param!(panel, 1, 3, 0.:0.005:0.5, 0.01,
            v -> "Modal damping : $v",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "0.1", "0.2", "0.3", "0.4", "0.5"], offset=2)

        f1 = add_param!(panel, 2, 1, 1.:0.1:20., 5.,
            v -> "Natural freq. 1 : $v Hz",
            [0., 9., 19.]/19., ["1", "10", "20"], offset=2)

        x01 = add_param!(panel, 2, 2, -1.:0.1:1., 0.,
            v -> "Initial displacement x₁ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"], offset=2)

        v01 = add_param!(panel, 2, 3, -10.:0.5:10., 0.,
            v -> "Initial velocity v₁ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"], offset=2)

        f2 = add_param!(panel, 3, 1, 1.:0.1:20., 5.,
            v -> "Natural freq. 2 : $v Hz",
            [0., 9., 19.]/19., ["1", "10", "20"], offset=2)

        x02 = add_param!(panel, 3, 2, -1.:0.1:1., 0.,
            v -> "Initial displacement x₂ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"], offset=2)

        v02 = add_param!(panel, 3, 3, -10.:0.5:10., 0.,
            v -> "Initial velocity v₂ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"], offset=2)

        # Masses and natural frequencies of the coupled system
        fn = lift((f1, f2) -> coupled_frequencies(m1, m2, f1, f2), f1, f2)
        Label(panel[1:2, 3], lift(fn -> "m₁ = $m1 kg, m₂ = $m2 kg\nModes : $(round(fn[1], digits = 2)) et $(round(fn[2], digits = 2)) Hz", fn),
            halign=:left, tellwidth=false, justification=:left, fontsize=15., color=:gray20)

        # Checkbox to freeze the y-axis scale
        option = GridLayout(panel[9, 1:3], halign=:left)
        checkbox = Checkbox(option[1, 1], checked=false, size=18., checkboxcolor_checked=colorant"#7fc4e8",
            checkboxstrokecolor_checked=colorant"#7fc4e8", checkboxstrokecolor_unchecked=:gray40)
        Label(option[1, 2], "Freeze y-axis", fontsize=16., color=:gray20)
        colgap!(option, 8.)
        y_fixed = checkbox.checked

        for i in 1:8
            rowgap!(panel, i, isodd(i) ? 8. : 30.)
        end
        colgap!(panel, 1, 24.)
        colgap!(panel, 2, 24.)

        colsize!(fig.layout, 1, Relative(0.5))
        colsize!(fig.layout, 2, Relative(0.5))

        # Forced response of the 2DOF system
        x = lift(response, f1, f2, ξ, A, k, f, x01, x02, v01, v02)
        x1 = lift(xs -> xs[1, :], x)
        x2 = lift(xs -> xs[2, :], x)

        hlines!(ax_x, 0., color=:gray30, linewidth=1.)
        lines!(ax_x, t, x1, linewidth=2., label="x₁(t)")
        lines!(ax_x, t, x2, linewidth=2., label="x₂(t)")
        xlims!(ax_x, t[begin], t[end])

        Legend(time_plot[1, 1], ax_x, orientation=:horizontal, halign=:left, framevisible=false, tellwidth=false, labelsize=16.)
        rowgap!(time_plot, 1, 4.)

        # Automatic scale, unless the scale is fixed (re-aligned on the curve at release)
        onany(x, y_fixed, update=true) do xs, fixed
            fixed && return
            xmin, xmax = extrema(xs)
            δ = max(0.05*(xmax - xmin), 1e-3)
            ylims!(ax_x, xmin - δ, xmax + δ)
        end

        # Component H_jk of the transfer function of the 2DOF system
        H = lift(transfer, f1, f2, ξ)
        Hjk = lift((Hs, j, k) -> Hs[j, k, :], H, j, k)
        amp = lift(Hs -> [isfinite(z) && z != 0 ? 20log10(abs(z)) : NaN for z in Hs], Hjk)
        phi = lift(Hs -> phase.(Hs), Hjk)

        name = lift((j, k) -> "H$(subscripts[j])$(subscripts[k])", j, k)
        connect!(ax_amp.title, lift((n, j, k) -> "Transfer function $n = x$(subscripts[j])/F$(subscripts[k])", name, j, k))
        connect!(ax_phi.title, lift(n -> "Phase of $n", name))

        # Locate the excitation frequency on the graphs
        f_exc = lift(f -> [f], f)

        for ax in (ax_amp, ax_phi)
            vlines!(ax, f_exc, color=colorant"#e36837", linewidth=1.5, linestyle=:dash)
            xlims!(ax, freq[begin], freq[end])
        end

        lines!(ax_amp, freq, amp, color=colorant"#1f77b4", linewidth=2.)
        lines!(ax_phi, freq, phi, color=colorant"#1f77b4", linewidth=2.)
        ylims!(ax_phi, -2π - 0.25, 0.25)

        onany(amp, y_fixed, update=true) do as, fixed
            fixed && return
            amin, amax = extrema(filter(isfinite, as))
            ylims!(ax_amp, 10floor(max(amin, amax - 150.)/10), 10ceil(amax/10))
        end

        fig
    end

    screen = display(fig)

    # Keep the window open when the script is run from a terminal
    if abspath(PROGRAM_FILE) == @__FILE__
        wait(screen)
    end
end

"""
    newmark_2dof(m1, m2, t)

Launch an interactive application to visualize the influence of Newmark parameters on the response of a two-degree-of-freedom (2DOF) system

**Inputs**
- `m1`: Mass of the first degree of freedom (default: 1.0)
- `m2`: Mass of the second degree of freedom (default: 1.0)
- `t`: Time vector over which the response is computed (default: 0.0 to 5.0 seconds with a step of 1e-3)

**Output**
- App
"""
function newmark_2dof(m1 = 1., m2 = 1., t = 0.:1e-3:5.)
    # Mass matrix
    M = mass(m1, m2)

    # Duration of the simulation
    T = t[end]

    # Special cases of the Newmark method : (button name, label on the stability diagram, γ, β)
    schemes = [("Central differences", "Central diff.", 0.5, 0.),
               ("Fox & Goodwin", "Fox & Goodwin", 0.5, 1/12),
               ("Linear acceleration", "Linear acc.", 0.5, 1/6),
               ("Constant average acceleration", "Average acc.", 0.5, 1/4)]

    # Matrice d'amortissement correspondant à l'amortissement modal forfaitaire
    function damping_matrix(K, ξ)
        ωn, Φn = eigenmode(K, M)
        return modal_damping_matrix(M, ωn, ξ, Φn)
    end

    force_vector(A, k) = k == 1 ? [A, 0.] : [0., A]

    # Reference solution
    function reference_response(p, ts)
        (f1, f2, ξ, A, k, f, x01, x02, v01, v02) = p
        prob = HarmonicModalTimeProblem(stiffness(m1, m2, f1, f2), M, ξ, force_vector(A, k), f, ([x01, x02], [v01, v02]), ts)
        return solve(prob).u
    end

    # Newmark solution and reference solution, computed on the same time grid
    function newmark_response(p, γ, β, Δt)
        (f1, f2, ξ, A, k, f, x01, x02, v01, v02) = p
        K = stiffness(m1, m2, f1, f2)
        ts = 0.:Δt:T
        F = force_vector(A, k)*cos.(2π*f*ts)'

        prob = DirectTimeProblem(K, M, damping_matrix(K, ξ), F, ([x01, x02], [v01, v02]), ts)

        return ts, solve(prob, Newmark(γ0=γ, β0=β), progress=false).u, reference_response(p, ts)
    end

    # Critical time step for the Newmark method
    function critical_step(f1, f2, ξ, γ, β)
        γ < 0.5 && return NaN
        β ≥ γ/2 && return Inf

        ωmax = 2π*maximum(coupled_frequencies(m1, m2, f1, f2))
        a = γ/2 - β
        b = ξ*(γ - 0.5)

        return (b + √(a + b^2))/(a*ωmax)
    end

    # Stability region of the Newmark method
    function stability_region(γ, β)
        γ < 0.5 && return "Unstable"
        β < γ/2 && return "Conditional stability"
        return β < (γ + 0.5)^2/4 ? "Unconditionally stable, real roots" : "Unconditionally stable, complex roots"
    end

    # Stability text for the Newmark method
    function stability_text(f1, f2, ξ, γ, β, Δt)
        region = stability_region(γ, β)
        Δtc = critical_step(f1, f2, ξ, γ, β)
        (isnan(Δtc) || isinf(Δtc)) && return region

        return "$region : Critical Δt = $(round(Δtc, sigdigits = 3)) s" * (Δt < Δtc ? ", Δt lower : stable" : ", Δt upper : divergence")
    end

    # Relative error between the Newmark solution and the reference solution
    function relative_error(x, xref)
        e = norm(x - xref)/max(norm(xref), 1e-12)
        return isfinite(e) && e < 10. ? "$(round(100e, sigdigits = 3)) %" : "> 1000 %"
    end

    # Format the points to be plotted (divergent values are not plotted)
    function curve_points(ts, xs, xmax)
        return Point2f[(ti, abs(xi) < 1e6*xmax ? xi : NaN) for (ti, xi) in zip(ts, xs)]
    end

    ## Application
    GLMakie.activate!(title="Influence of Newmark parameters on the response of a 2DOF system")
    theme = theme_choice(:makie, fonts=(regular="Lato", bold="Lato Bold"))

    fig = with_theme(theme) do
        fig = Figure(size=(1500, 900), figure_padding=30.)

        Label(fig[1, 1:2], "Influence of Newmark parameters on the response of a 2DOF system", fontsize=34., halign=:left, tellwidth=false, color=:gray20)

        axis_style = (xlabelsize=18., ylabelsize=18., xticklabelsize=14., yticklabelsize=14.,
            topspinevisible=false, rightspinevisible=false)

        # Graphiques
        time_plot = GridLayout(fig[2, 1])
        ax_x = Axis(time_plot[2, 1]; xlabel="Time (s)", ylabel="Displacement (m)", yticklabelspace=55., axis_style...)
        ax_s = Axis(fig[2, 2]; xlabel="γ", ylabel="β", xticks=0.:0.5:2.5, yticks=0.:0.5:1.5, axis_style...)

        # Contrôle panel
        frame = GridLayout(fig[3, 1:2])
        Box(frame[1, 1], color=:transparent, strokewidth=1., strokecolor=:gray30, cornerradius=6.)
        panel = GridLayout(frame[1, 1], alignmode=Outside(18., 18., 14., 14.))

        s_γ = add_slider!(panel, 1, 1, 0.:0.01:2.5, 0.5,
            v -> "γ : $v",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "0.5", "1", "1.5", "2", "2.5"])

        # The step 1/120 allows to reach exactly 1/12, 1/6 and 1/4
        s_β = add_slider!(panel, 2, 1, (0:180)/120, 0.25,
            v -> "β : $(round(v, digits = 3))",
            [0., 1/3, 2/3, 1.], ["0", "0.5", "1", "1.5"])

        s_Δt = add_slider!(panel, 3, 1, 0.001:0.001:0.2, 0.02,
            v -> "Time step Δt : $v s",
            [0., 0.25, 0.5, 0.75, 1.], ["0", "0.05", "0.1", "0.15", "0.2"])

        f1 = add_param!(panel, 1, 2, 1.:0.5:20., 5.,
            v -> "Natural freq. of syst. 1 : $v Hz",
            [0., 9., 19.]/19., ["1", "10", "20"])

        f2 = add_param!(panel, 2, 2, 1.:0.5:20., 5.,
            v -> "Natural freq. of syst. 2 : $v Hz",
            [0., 9., 19.]/19., ["1", "10", "20"])

        ξ = add_param!(panel, 3, 2, 0.:0.005:0.5, 0.01,
            v -> "Modal damping : $v",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "0.1", "0.2", "0.3", "0.4", "0.5"])

        A = add_param!(panel, 1, 3, 0.:10.:500., 100.,
            v -> "Amplitude F : $(round(Int, v)) N",
            [0., 0.2, 0.4, 0.6, 0.8, 1.], ["0", "100", "200", "300", "400", "500"])

        f = add_param!(panel, 2, 3, 0.5:0.5:30., 3.,
            v -> "Excitation freq. : $v Hz",
            [0., 9.5, 19.5, 29.5]/29.5, ["0.5", "10", "20", "30"])

        Label(panel[5, 3], "Force applied on", halign=:left, tellwidth=false, fontsize=16., color=:gray20)
        menu_k = Menu(panel[6, 3], options=[("Mass 1", 1), ("Mass 2", 2)], default="Mass 2", fontsize=15., direction=:up)
        k = menu_k.selection

        x01 = add_param!(panel, 1, 4, -1.:0.1:1., 0.,
            v -> "Initial displacement x₁ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"])

        v01 = add_param!(panel, 2, 4, -10.:0.5:10., 0.,
            v -> "Initial velocity v₁ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"])

        x02 = add_param!(panel, 1, 5, -1.:0.1:1., 0.,
            v -> "Initial displacement x₂ : $v m",
            [0., 0.5, 1.], ["-1", "0", "1"])

        v02 = add_param!(panel, 2, 5, -10.:0.5:10., 0.,
            v -> "Initial velocity v₂ : $v m/s",
            [0., 0.25, 0.5, 0.75, 1.], ["-10", "-5", "0", "5", "10"])

        # Special cases
        presets = GridLayout(panel[7, 1:5])
        Label(presets[1, 1], "Special cases :", fontsize=16., color=:gray20)
        for (i, (name, _, γp, βp)) in enumerate(schemes)
            button = Button(presets[1, i+1], label=name, fontsize=15., buttoncolor=colorant"#e9e9e9")
            on(button.clicks) do _
                set_close_to!(s_γ, γp)
                set_close_to!(s_β, βp)
            end
        end
        Box(presets[1, 6], visible=false)

        # Case à cocher de blocage de l'échelle des ordonnées
        option = GridLayout(presets[1, 7], halign=:left)
        checkbox = Checkbox(option[1, 1], checked=false, size=18., checkboxcolor_checked=colorant"#7fc4e8",
            checkboxstrokecolor_checked=colorant"#7fc4e8", checkboxstrokecolor_unchecked=:gray40)
        Label(option[1, 2], "Freeze y-axis", fontsize=16., color=:gray20)
        colgap!(option, 8.)
        y_fixed = checkbox.checked
        colgap!(presets, 12.)

        γ, β, Δt = s_γ.value, s_β.value, s_Δt.value

        for i in 1:5
            rowgap!(panel, i, isodd(i) ? 6. : 28.)
        end
        rowgap!(panel, 6, 30.)
        for i in 1:4
            colgap!(panel, i, 26.)
        end

        colsize!(fig.layout, 1, Relative(0.5))
        colsize!(fig.layout, 2, Relative(0.5))

        # Numerical and reference solutions of the 2DOF system
        p = lift(tuple, f1, f2, ξ, A, k, f, x01, x02, v01, v02)
        sol = lift(newmark_response, p, γ, β, Δt)
        xmax = lift(s -> max(maximum(abs, s[3]), 1e-3), sol)

        hlines!(ax_x, 0., color=:gray30, linewidth=1.)

        for (j, name, color) in ((1, "x₁", Cycled(1)), (2, "x₂", Cycled(2)))
            lines!(ax_x, lift((s, xm) -> curve_points(s[1], s[2][j, :], xm), sol, xmax), color=color, linewidth=2., linestyle=:solid, label="$name Newmark")
            lines!(ax_x, lift((s, xm) -> curve_points(s[1], s[3][j, :], xm), sol, xmax), color=color, alpha=0.6, linewidth=1.5, linestyle=:dash, label="$name reference")
        end
        xlims!(ax_x, 0., T)

        Legend(time_plot[1, 1], ax_x, orientation=:horizontal, halign=:left, framevisible=false, tellwidth=false, labelsize=16.)
        rowgap!(time_plot, 1, 4.)

        # The axis follows the numerical response up to 5 times the analytical amplitude, unless the scale is fixed
        onany(sol, xmax, y_fixed) do s, xm, fixed
            fixed && return
            finite = filter(isfinite, s[2])
            ylim = 1.05*min(isempty(finite) ? Inf : maximum(abs, finite), 5xm)
            ylims!(ax_x, -max(ylim, 1.05xm), max(ylim, 1.05xm))
        end
        notify(sol)

        # Stability diagram (βmin leaves visible the point of centered differences)
        βmin = -0.05
        γs = range(0.5, 2.5, length=200)
        band!(ax_s, [0., 0.5], [βmin, βmin], [1.5, 1.5], color=colorant"#e4e4e4")
        band!(ax_s, γs, fill(βmin, length(γs)), γs/2, color=colorant"#fde3e4")
        band!(ax_s, γs, γs/2, min.((γs .+ 0.5) .^ 2/4, 1.5), color=colorant"#e4e4fb")

        γl = range(0., 2.5, length=200)
        lines!(ax_s, γl, γl/2, color=:red, linewidth=3., linestyle=:solid)
        lines!(ax_s, γl, (γl .+ 0.5) .^ 2/4, color=:blue, linewidth=3., linestyle=:solid)
        vlines!(ax_s, 0.5, color=:black, linewidth=1.5, linestyle=:solid)

        region_style = (align=(:center, :center), fontsize=15., color=:gray20)
        text!(ax_s, 0.25, 0.75; text="Instable", region_style...)
        text!(ax_s, 1.0, 1.2; text="Unconditional stability\nComplex roots", region_style...)
        text!(ax_s, 2.18, 1.35; text="Unconditional\nstability\nReal roots", region_style...)
        text!(ax_s, 2.05, 0.3; text="Conditional stability", region_style...)
        text!(ax_s, 1.13, 0.75; text="β = (γ + ½)²/4", align=(:center, :center), fontsize=15., color=:blue, rotation=0.72)
        text!(ax_s, 1.72, 0.78; text="β = γ/2", align=(:center, :center), fontsize=15., color=:red, rotation=0.46)

        # Special cases and current point
        scatter!(ax_s, [Point2f(s[3], s[4]) for s in schemes], color=:black, markersize=9.)
        for (_, label, γp, βp) in schemes
            text!(ax_s, γp + 0.04, βp, text=label, align=(:left, :center), fontsize=13., color=:gray20)
        end
        scatter!(ax_s, lift(Point2f, γ, β), color=:red, markersize=16., strokecolor=:black, strokewidth=1.5)

        limits!(ax_s, 0., 2.5, βmin, 1.5)
        connect!(ax_s.title, lift((γ, β) -> "(γ, β) = ($γ, $(round(β, digits = 3)))", γ, β))

        # Summary: stability, deviation from the analytical solution and system data
        error_text = lift(s -> relative_error(s[2], s[3]), sol)
        status = lift((f1, f2, ξ, γ, β, Δt, e) -> stability_text(f1, f2, ξ, γ, β, Δt) * "\nDeviation from the analytical solution : $e", f1, f2, ξ, γ, β, Δt, error_text)
        system = lift(fn -> "Modes : $(round(fn[1], digits = 2)) et $(round(fn[2], digits = 2)) Hz, m₁ = $m1 kg, m₂ = $m2 kg", lift((f1, f2) -> coupled_frequencies(m1, m2, f1, f2), f1, f2))

        info = GridLayout(panel[5:6, 4:5])
        Label(info[1, 1], status, halign=:left, tellwidth=false, justification=:left, fontsize=15., color=:gray20)
        Label(info[2, 1], system, halign=:left, tellwidth=false, fontsize=15., color=:gray20)
        rowgap!(info, 1, 2.)

        fig
    end

    screen = display(fig)

    # Keep the window open when the script is run from a terminal
    if abspath(PROGRAM_FILE) == @__FILE__
        wait(screen)
    end
end
