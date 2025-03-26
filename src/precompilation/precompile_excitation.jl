# Excitation
rectangle = Rectangle(F₀, tstart, duration)
triangle = Triangle(F₀, tstart, duration)
hammer = Hammer(F₀, tstart, 8.7, 6e-4)
smoothrect = SmoothRect(F₀, tstart, duration, 5Δt)
sinewave = SineWave(F₀, tstart, duration, 10.)
halfsine = HalfSine(F₀, tstart, duration)
haversine = HaverSine(F₀, tstart, duration)
sweptsine = SweptSine(F₀, tstart, duration, 10., 20.)
gaussianpulse = GaussianPulse(F₀, tstart, duration, 10.)
colorednoise = ColoredNoise(F₀, tstart, duration, color = :pink)

@compile_workload begin
    Ft = excitation(rectangle, t)
    Ft .= excitation(triangle, t)
    Ft .= excitation(hammer, t)
    Ft .= excitation(smoothrect, t)
    Ft .= excitation(sinewave, t)
    Ft .= excitation(halfsine, t)
    Ft .= excitation(haversine, t)
    Ft .= excitation(sweptsine, t)
    Ft .= excitation(gaussianpulse, t)
    Ft .= excitation(colorednoise, t)
end