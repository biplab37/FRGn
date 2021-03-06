#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. This is to reproduce the result by Bauer et al.                        ###
#########################################################################################

using FRGn

## Initialisation

const m = 301 # number of cutoffs
const n = 217 # number of momenta

velocity = zeros(n,m)
dielectric = zeros(n,m)

@doc raw"""
    velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)

This function returns the integrand of the FRG equation for the velocity renormlisation.
## Args
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
    ## Theta function implementation with conditional

    k = sqrt(cutoff^2 + momentum^2 - 2*cutoff*momentum*cos(2.0*phi))

    epsilon = get_dielectric(dielectric,k,cutoff,m,n)

    if k==0
        return 0.0
    else
        return 2.2*cos(2.0*phi)*cutoff/(pi*epsilon*momentum*k)
    end
end

@doc raw"""
    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)

This function returns the integrand of the FRG equation for the dielectric function renormlisation.
## Args
    velocity   (Array) : The velocity array
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, cutoff, m, n)

        return 4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)*sqrt((k1+k2)^2 - momentum^2))
    end
end

## Boundary values initialisation
velocity[:,m] .= 1.0
dielectric[:,m] .= 1.0

## solving exact FRG using an user defined function from FRGn Package

rg_procedure(velocity,dielectric,velocity_integrand, dielectric_integrand ,m,n)


## Plots using user defined functions in FRGn Package

plot_velocity(velocity[:,1],"bauer_renormalised_velocity.pdf")

plot_dielectric(dielectric[:,1],"bauer_renormalised_dielectric.pdf")

## Save the data for future usage
# using JLD
# save("bauer.jld","velocity",velocity,"dielectric",dielectric)

# using HDF5
# h5write("bauer.h5","velocity",velocity)
# h5write("bauer.h5","dielectric",dielectric)
