# NLControl
***Under development***

A Julia package to solve optimal control problems for Non-Linear systems of differential equations.

[![Build Status](https://github.com/qmarcou/NLControl.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/qmarcou/NLControl.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/qmarcou/NLControl.jl/graph/badge.svg?token=DGOAEE7HNP)](https://codecov.io/gh/qmarcou/NLControl.jl)

## Aim
This package aims to enable the use of JuMP.jl solvers to solve optimal control problems while using a simple API akin to DifferentialEquations.jl to encode the dynamical system information. On top of a simplified interface to solve optimal control problems, NLControl.jl aim to provide access to higher order Runge-Kutta integration schemes than Euler integration and solution checking using DifferentialEquations.jl adaptive integration.  

## Implementation TODO
- automated solution checking
- RK integration
- handling of multiple controls

## Limitations
- Optimization of control time steps is not possible with the current implmentation.