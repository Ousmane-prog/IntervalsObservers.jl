# IntervalsObservers

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ousmane-prog.github.io/IntervalsObservers.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ousmane-prog.github.io/IntervalsObservers.jl/dev/)
[![Build Status](https://github.com/ousmane-prog/IntervalsObservers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ousmane-prog/IntervalsObservers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ousmane-prog/IntervalsObservers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ousmane-prog/IntervalsObservers.jl)

# IntervalObservers.jl

A Julia package for designing and simulating interval-based observers for linear and nonlinear systems with model uncertainty.

## Overview

`IntervalObservers.jl` provides tools for:

- **Linear Observer Design**: Compute observer gains for linear systems with guaranteed convergence
- **Nonlinear Systems**: Handle systems with interval-bounded nonlinear uncertainty terms
- **State Estimation**: Obtain upper and lower bounds on system states despite model uncertainty
- **Visualization**: Plot state interval estimates and trajectories

## Key Features

- Automatic gain computation using pole placement
- Support for Metzler matrices and monotone dynamics
- Michaelis-Menten type nonlinearities
- Change of basis transformations for system stabilization
- Comprehensive error checking and validation

## Installation

```julia
using Pkg
Pkg.add("IntervalObservers")
```

## Quick Start

See the [Nonlinear System Example](nonlinear_system_example.md) for a complete walkthrough.

## Documentation

This documentation includes:

- **[Nonlinear System Example](nonlinear_system_example.md)** - Complete example showing how to design and simulate an interval observer for a nonlinear system

