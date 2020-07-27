#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

using FRGn
using QuadGK

## Initialisation

const m = 329 # number of cutoffs
const n = 343 # number of momenta

velocity = zeros(n,m)
dielectric = zeros(n,m)
dielectric2 = zeros(n,m)

@doc raw"""
    velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)

This function returns the integrand of the FRG equation for the velocity renormlisation.
## Args
    velocity   (Array) : The velocity array
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate (should run from 0 to pi/2)
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
    ## Theta function implementation with conditional
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        eps1,eps2 = get_dielectric(dielectric, k1, k2, m, n, i)

        return 2.2*((momentum^2 - k1^2 + k2^2)/(momentum^2*eps1) + (momentum^2 + k1^2 - k2^2)/(momentum^2*eps2))/(2.0*pi*sqrt((k1+k2)^2 - momentum^2))
    end
end

@doc raw"""
    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)

This function returns the integrand of the FRG equation for the dielectric function renormlisation.
## Args
    velocity   (Array) : The velocity array
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
"""
function dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, m, n, i)

        return 4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)*sqrt((k1+k2)^2 - momentum^2))
    end
end

function dielectric2_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, m, n, i)

        return -4.4*momentum*sin(phi)^2/(pi*(k1*vel1 + k2*vel2)^3*sqrt((k1+k2)^2 - momentum^2))
    end
end

function rg_procedure2(velocity::Array{Float64,2}, dielectric::Array{Float64,2},dielectric2::Array{Float64,2}, velocity_integrand, dielectric_integrand,dielectric2_integrand, m::Int64, n::Int64)

    dcutoff::Float64 = 1.0/m

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi,m,n,i)
            dielectric_integrand_phi(phi) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi,m,n,i)
            dielectric2_integrand_phi(phi) = dielectric2_integrand(velocity,dielectric,momentum,cutoff,phi,m,n,i)
            ## Solving ODE's using Euler method
            velocity[j,m-i+1] = velocity[j,m-i+2] + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            dielectric[j,m-i+1] = dielectric[j,m-i+2] + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            dielectric2[j,m-i+1] = dielectric2[j,m-i+2] + dcutoff*quadgk(dielectric2_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
        end
    end
end

## Boundary values initialisation
velocity[:,m] .= 1.0
dielectric[:,m] .= 1.0
dielectric2[:,m] .= 0.0

## solving exact FRG using an user defined function from FRGn Package

rg_procedure2(velocity,dielectric,dielectric2,velocity_integrand, dielectric_integrand,dielectric2_integrand ,m,n)


## Plots using user defined functions in FRGn Package

plot_velocity(velocity[:,1])

plot_dielectric(dielectric[:,1])
plot_dielectric(dielectric2[:,1])

## Save the data for future usage
using JLD
save("correction.jld","velocity",velocity,"dielectric",dielectric,"dielectric2",dielectric2)

# using HDF5
# h5write("coupled.h5","velocity",velocity)
# h5write("coupled.h5","dielectric",dielectric)
