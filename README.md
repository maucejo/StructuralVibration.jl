# StructuralVibration.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://maucejo.github.io/StructuralVibration.jl/)
[![Generic badge](https://img.shields.io/badge/Version-0.1.0-cornflowerblue.svg)]()
[![MIT License](https://img.shields.io/badge/License-MIT-forestgreen)](https://github.com/maucejo/elsearticle/blob/main/LICENSE)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This Julia package is intended to provide a set of tools for generating vibration data either to validate new numerical methods or to serve as input data for solving inverse problems such as input state estimation via Bayesian filtering, Bayesian regularization, machine learning, etc.

The package is under active development and is not yet registered in the Julia General Registry. The package is developed as part of my research activities on source identification in structural dynamics and acoustics.

## Installation

The package is not yet registered in the Julia General Registry. To install the package, you can use the following command:

```julia
(Yourenv) pkg> add "git@https://github.com/maucejo/StructuralVibration.jl.git"
```

## Package features

The package is organized into the following functionalities:

### Mechanical models

- **Discrete models**
    - [x] Spring-mass-damper SDOF system
    - [ ] Spring-mass-damper MDOF system
    - [ ] FE models of bar, rod and beams
- **Modal models of continuous model**
    - [x] Longitudinal and torsional bars for various boundary conditions
        - Clamped-clamped
        - Clamped-free
        - Free-free
    - [x] Euler-Bernoulli beams for various boundary conditions
        - Simply-supported at both ends
        - Clamped-clamped
        - Clamped-simply-supported
        - Clamped-free
        - Simply-supported-free
    - [x] Simply-supported plate

### Vibration data generation

- **Excitation signals**
    - [x] White noise
    - [x] Triangle wave
    - [x] Rectangular wave
    - [x] Hammer impact
    - [x] Smoothed rectangular wave
    - [ ] Chirp

- **Solution for SDOF systems**
    - [x] Free response
    - [x] Forced response due to a harmonic force or a base motion
    - [x] Forced response due to any external force

- **Time-domain integration schemes for linear systems**
    - Central difference scheme
    - RK4
    - Newmark-beta method
    - HHT
    - WBZ
    - Generalized-alpha
    - Mid-Point rule

- **Frequency-domain calculations for linear systems**
    - Frequency spectrum
        - [x] Modal summation
        - [ ] Direct method
    - Frequency response function (FRF)
        - [x] Modal summation
        - [ ] Direct method

### Measurement noise

- [x] Addition of Gaussian white noise with a prescribed SNR
- [x] Measurement noise variance estimation algorithms from noisy data
    - [x] Bayesian estimation
    - [x] D'Errico's method  - [Link to the Matlab version](https://fr.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise)
- [x] SNR estimation from estimated measurement noise variance