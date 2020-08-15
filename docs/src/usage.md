# User Guide

```@contents
Pages = ["usage.md"]
```

To be able  solve exact FRG equations one need to follow the steps below.

## Preparation

Before starting coding one need to do the necessary analytical calculations in order to obtain the exact FRG equations that one wants to solve. Then one need to simplify the equations such as to mould them into the following form

$\partial_{\Lambda}f_i = \int_0^{\pi/2}\hat{f_i}$

where $f_i$'s are the observables that we want to solve and $\hat{f_i}$ are the integrands corrosponding to the observables. Note the limit of the integration is from $0$ to $\pi/2$. If one needs to solve a problem which can not be turned into this form please raise an issue. Once we have the closed set of equations we are ready to solve them.

## Initialisation

Note that the observables $f_i$'s are dependent on the momentum and the running cutoff. So we would like them to store in a 2-dimensional array. So the code would look something like

```julia
using FRGn

m = 123 # no of cutoffs
n = 456 # no of momenta
observable1 = zeros(n,m)
observable2 = zeros(n,m)
...

```


## Define integrands

Once the observables are defined the integrands on the rhs of FRG equation have to be provided. This should be done as a function of the following form

```julia
function integrand1(observable1 ,observable2, momentum, cutoff, phi, m, n) 
	...
end
```

## Solving

Before solving the differential equations we need to supply the boundary condition. Boundary condition are comes from the fact that we know the observables at the highest cutoff i.e. before starting the renormalisation group procedure. For example

```julia
observable1[:,m] .= 1
```
where we redefined the variable so that the value of observable was unity at the highest cutoff.

Finally, to solve the FRG equations we call the *rg_procedure* funciton

```julia
rg_procedure(observable1, observable2, integrand1, integrand2, m, n)
```

This will update the observables in place.

## Visualise Results

To visulaise the result one can plot the resultant observables for specific cases using *Plots* or *PyPlot*. Here we provide two simple function to plot the renormlised velocity and dielectric function as follows
```julia
plot_velocity(velocity[:,1],"renormalised_velocity.png")
``` 
this will also save the figure in the given name. Omiting that variable will not save the figure.