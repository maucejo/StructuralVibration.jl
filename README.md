# StructuralVibration.jl

[![Generic badge](https://img.shields.io/badge/Version-1.0.0-cornflowerblue.svg)]()
[![MIT License](https://img.shields.io/badge/License-MIT-forestgreen)](https://github.com/maucejo/elsearticle/blob/main/LICENSE)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DispatchDoctor](https://img.shields.io/badge/%F0%9F%A9%BA_tested_with-DispatchDoctor.jl-blue?labelColor=white)](https://github.com/MilesCranmer/DispatchDoctor.jl)

`StructuralVibration.jl` is a Julia package designed to generate vibration data for mechanical systems. This package can be used for educational and research purposes. As a teacher, you can use it to illustrate vibration theory and its practical implementation. As a researcher, you can use it to test new methods on synthetic data.

::: {.callout-warning}
The package is under active development and is developed as part of my research activities on source identification in structural dynamics and acoustics.
:::

## Features

The package provides the following features:

### Mechanical models

- **Discrete models**
    - Spring-mass-damper SDOF system
    - Spring-mass-damper MDOF system
    - FE models of bar, rod, strings and beams
- **Continuous models**
    - Longitudinal bars
    - Torsional bars
    - Strings
    - Beams
    - Rectangular plates
    - Rectangular membranes
- **State space model**
    - Continuous state-space representation
    - Discrete state-space representation
        - Zero-order hold (ZOH)
        - First-order hold (FOH)
        - Band-limited hold (BLH)
        - RK4

### Vibration data generation

- **Excitation signals**
    - Rectangular wave
    - Triangular wave
    - Hammer impact
    - Smoothed rectangular wave
    - Sine wave
    - Half-sine pulse
    - Harversine pulse
    - Swept sine wave
    - Gaussian pulse
    - Colored noise

- **Solution for SDOF systems**
    - Free response
    - Forced response due to a harmonic force or a base motion
    - Forced response due to any external force or base motion (Duhamel's integral)

- **Time-domain integration schemes for linear second order systems**
    - Central difference scheme
    - RK4
    - Newmark-beta method
    - Linear acceleration method
    - Fox-Goodwin method
    - HHT
    - WBZ
    - Generalized-alpha
    - Mid-Point rule

- **Frequency-domain calculations for linear systems**
    - Frequency spectrum
        - Modal summation
        - Direct method
    - Frequency response function (FRF)
        - Modal summation
        - Direct method

- **State-space solvers**
    - Time domain
        - RK4 for continuous systems
        - ZOH, FOH, BLH, RK4 for discrete models
    - Frequency spectrum
        - Modal summation
        - Direct method
    - Frequency response function (FRF)
        - Modal summation
        - Direct method

- **Measurement noise**
    - Additive Gaussian white noise (agwn) with a prescribed SNR
    - Additive Colored noise (acn)
    - Multiplicative noise
    - Mixed - agwn + multiplicative noise

### Signal processing

- **Measurement noise variance estimation algorithms from noisy data**
    - Regularization-based estimation
    - D'Errico's method  - [Link to the Matlab version](https://fr.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise)
    - SNR estimation from estimated measurement noise variance

- **Denoising algorithms**
    - Regularization
    - Kalman filtering

- **Signal preprocessing**
    - Detrending data using polynomial fit
    - Gradient calculation using interpolation

- **Signal estimation**
    - Transfer functions estimation (H1, H2, H3, Hv)
    - Welch method (PSD, ESD, Autopower, Autopower linear)
    - Signal spectrum estimation

### Modal analysis

- **Modal extraction - SDOF methods**
    - Peak picking method
    - Circle fit method
    - Least-squares fit method

- **Modal extraction - MDOF methods**
    - Least-squares complex exponential (LSCE)
    - Least-squares complex frequency (LSCF)
    - Polyreference least-squares complex frequency (PLSCF)
    - Stabilization diagram

- **Analysis indicators**
    - Modal Overlap Factor (MOF)
    - Mode overcomplexity value (MOV)
    - Mode Phase Collinearity (MPC)
    - Mode Complexity Factor (MCF)
    - Mean Phase Deviation (MPD)
    - Modal Assurance Criterion (MAC)
    - Coordinate Modal Assurance Criterion (COMAC)
    - Enhanced Coordinate Modal Assurance Criterion (ECOMAC)
    - Frequency Response Assurance Criterion (FRAC)
    - Complex Mode Indicator Function (CMIF)
    - Power Spectrum indicator function (PSIF)

- **Modal synthesis**
    - Estimation of the FRF residuals
    - Reconstruct FRF from modal parameters

### Visualization

- Bode plot
- 2D and 3D Nyquist plot
- Waterfall plot
- General 2D plot
