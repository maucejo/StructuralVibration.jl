"""
    add_marks!(slider, fracs, labels; fontsize = 13., color = :gray40)

Display of graduations under a slider

**Inputs**
- `slider`: Slider object
- `fracs`: Fractions of the slider range where the graduations should be placed (between 0 and 1)
- `labels`: Labels corresponding to the graduations
- fontsize: Font size of the labels (default: 13)
- color: Color of the labels (default: :gray40)

**Output**
- Text objects corresponding to the graduations
"""
function add_marks!(slider, fracs, labels; fontsize = 13., color = :gray40)
    positions = lift(slider.layoutobservables.computedbbox) do bb
        h = height(bb)
        Point2f[(left(bb) + h/2 + f*(width(bb) - h), bottom(bb) - 5.) for f in fracs]
    end
    text!(slider.blockscene, positions, text = labels, align = (:center, :top), fontsize = fontsize, color = color)
end

"""
    add_slider!(panel, row, col, range, startvalue, format, fracs, labels; deps = (), offset = 0)

Add a parameter to the control panel (current value + slider + graduations) and return its slider

Use it instead of `add_param!` when the slider itself is needed, e.g. to set its value with `set_close_to!`

**Inputs**
- See `add_param!`

**Output**
- `slider`: Slider object of the parameter
"""
function add_slider!(panel, row, col, range, startvalue, format, fracs, labels; deps = (), offset = 0)
    label_row = 2row - 1 + offset

    slider = Slider(panel[label_row + 1, col], range = range, startvalue = startvalue, linewidth = 12.,
                    color_active = colorant"#7fc4e8", color_active_dimmed = colorant"#b5dcf1", color_inactive = colorant"#e9e9e9")
    Label(panel[label_row, col], lift(format, slider.value, deps...), halign = :left, fontsize = 16., color = :gray20)
    add_marks!(slider, fracs, labels)

    return slider
end

"""
    add_param!(panel, row, col, range, startvalue, format, fracs, labels; deps = (), offset = 0)
    add_param!(panel, row, range, startvalue, format, fracs, labels; deps = (), offset = 0)

Add a parameter to the control panel: current value + slider + graduations

**Inputs**
- `panel`: GridLayout object corresponding to the control panel
- `row`: Index of the parameter row in the panel (each parameter row spans two rows of the grid)
- `col`: Column of the panel where the parameter should be added (default: 1)
- `range`: Range of the parameter (as a UnitRange or StepRange)
- `startvalue`: Initial value of the parameter
- `format`: Function to format the current value of the parameter for display
- `fracs`: Fractions of the slider range where the graduations should be placed (between 0 and 1)
- `labels`: Labels corresponding to the graduations
- `deps`: Additional observables passed to `format` after the value of the parameter (default: `()`)
- `offset`: Number of grid rows reserved above the parameters (default: 0)

**Output**
- `slider.value`: Observable corresponding to the current value of the parameter
"""
add_param!(panel, row, col, range, startvalue, format, fracs, labels; kwargs...) = add_slider!(panel, row, col, range, startvalue, format, fracs, labels; kwargs...).value

add_param!(panel, row, range, startvalue, format, fracs, labels; kwargs...) = add_param!(panel, row, 1, range, startvalue, format, fracs, labels; kwargs...)

"""
    stiffness(m1, m2, f1, f2)
Compute the stiffness matrix of a 2-DOF system given its natural frequencies

**Inputs**
- `m1`: Mass of the first degree of freedom
- `m2`: Mass of the second degree of freedom
- `f1`: Natural frequency of the first mode (in Hz)
- `f2`: Natural frequency of the second mode (in Hz)

**Output**
- Stiffness matrix of the 2-DOF system
"""
function stiffness(m1, m2, f1, f2)
    k1 = m1*(2π*f1)^2
    k2 = m2*(2π*f2)^2

    return [k1 + k2 -k2; -k2 k2]
end

"""
    mass(m1, m2)
Compute the mass matrix of a 2-DOF system given its masses

**Inputs**
- `m1`: Mass of the first degree of freedom
- `m2`: Mass of the second degree of freedom

**Output**
- Mass matrix of the 2-DOF system
"""
function mass(m1, m2)
    return [m1 0.; 0. m2]
end

"""
    coupled_frequencies(f1, f2)

Compute the coupled natural frequencies of a 2-DOF system given its uncoupled natural frequencies

**Inputs**
- `m1`: Mass of the first degree of freedom
- `m2`: Mass of the second degree of freedom
- `f1`: Uncoupled natural frequency of the first mode (in Hz)
- `f2`: Uncoupled natural frequency of the second mode (in Hz)

**Output**
- Coupled natural frequencies of the 2-DOF system (in Hz)
"""
coupled_frequencies(m1, m2, f1, f2) = eigenmode(stiffness(m1, m2, f1, f2), mass(m1, m2))[1]/2π
