# User Guide

This package provides a bunch of functions listed below

```@index
Order = [:function]
```
## Initialisation


## Define integrands

To use these functions one have to first analytically find the exact FRG equation which they want to numerically solve. Usually the FRG equaitons are coupled integro-differential equation. To use the equations here one need to cast them in the following form

$\partial_{\Lambda}f_i = \int_0^{\pi/2}\hat{f_i}$

where $f_i$'s are the functions that we want to solve and $\hat{f_i}$ are the integrands corrosponding to the functions. Note the limit of the integration is from $0$ to $\pi$. The user need to provide the integrands as a function. The funciton should have following form

```julia
    integrand(functions::Array{Array{Float64,2},1}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
```
Note the following from the above function

 - The functions on the lhs of the FRG equations should first be initialised as a two dimensional array of size m $\times$ n. This is because we want to solve the equations in a grid with m cutoffs and n momenta (The cutoff and momenta is normalised with the highest cutoff $\Lambda_0$).
 - _phi_ is the angular coordinate that is to be integrated out.

## Solving

## Visualise Results