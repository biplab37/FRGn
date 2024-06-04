#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

using FRGn

## Initialisation

const m = 329 # number of cutoffs
const n = 343 # number of momenta

velocity = zeros(n, m)
dielectric = zeros(n, m)
dielectric2 = zeros(n, m)

@doc raw"""
    velocity_integrand(functions, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)

This function returns the integrand of the FRG equation for the velocity renormlisation.
## Args
    functions  (Array) : An Array of all the functions whose derivatives are to be solved.
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate (should run from 0 to pi/2)
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
"""
function velocity_integrand(functions, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
    ## Theta function implementation with conditional
    if cos(phi) <= 1 - 2 * cutoff / momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi) * momentum

        eps1, eps2 = get_dielectric(functions[2], k1, k2, cutoff, m, n)

        return 2.2 * ((momentum^2 - k1^2 + k2^2) / (momentum^2 * eps1) + (momentum^2 + k1^2 - k2^2) / (momentum^2 * eps2)) / (2.0 * pi * sqrt((k1 + k2)^2 - momentum^2))
    end
end

@doc raw"""
    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)

This function returns the integrand of the FRG equation for the dielectric function renormlisation.
## Args
    functions  (Array) : An Array of all the functions whose derivatives are to be solved.
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
"""
function dielectric_integrand(functions, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
    ## Theta function implementation
    if cos(phi) <= 1 - 2 * cutoff / momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi) * momentum

        vel1, vel2 = get_velocity(functions[1], k1, k2, cutoff, m, n)

        return 4.4 * momentum * sin(phi)^2 / (pi * (k1 * vel1 + k2 * vel2) * sqrt((k1 + k2)^2 - momentum^2))
    end
end

@doc raw"""
    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)

This function returns the integrand of the FRG equation for the correction to the dielectric function (order ω^2) renormlisation.
## Args
    functions  (Array) : An Array of all the functions whose derivatives are to be solved.
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
"""
function dielectric2_integrand(functions, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)
    if cos(phi) <= 1 - 2 * cutoff / momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi) * momentum

        vel1, vel2 = get_velocity(functions[1], k1, k2, cutoff, m, n)

        return -4.4 * momentum * sin(phi)^2 / (pi * (k1 * vel1 + k2 * vel2)^3 * sqrt((k1 + k2)^2 - momentum^2))
    end
end

## Boundary values initialisation
velocity[:, m] .= 1.0
dielectric[:, m] .= 1.0
dielectric2[:, m] .= 0.0

## solving exact FRG using an user defined function from FRGn Package

rg_procedure([velocity, dielectric, dielectric2], [velocity_integrand, dielectric_integrand, dielectric2_integrand], m, n)


plot(range(0, 1, n), -dielectric[:, 1] ./ dielectric2[:, 1], lab="", title="Plasmon Frequency", xlabel=L"k/\Lambda_0", ylabel=L"\omega_p/\Lambda_0", l=2)
savefig("work/plots/plasmon_frequency.pdf")
## Plots using user defined functions in FRGn Package

plot_velocity(velocity[:, 1])

plot_dielectric(dielectric[:, 1])
plot_dielectric(dielectric2[:, 1])

## Save the data for future usage
# using JLD
# save("correction.jld","velocity",velocity,"dielectric",dielectric,"dielectric2",dielectric2)

# using HDF5
# h5write("coupled.h5","velocity",velocity)
# h5write("coupled.h5","dielectric",dielectric)
