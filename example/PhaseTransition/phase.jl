#########################################################################################
###    This code solves the exact FRG equations numerically for graphene near Dirac   ###
###    points. Here we have introduced a cutoff for the Bososnic momenta as well.     ###
#########################################################################################

using FRGn

## Initialisation

const m = 239 # number of cutoffs
const n = 233 # number of momenta

println("Number of cutoffs taken: $m")
println("Number of momenta taken: $n")

velocity = zeros(n,m)
dielectric = zeros(n,m)

function velocity_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, temp::Float64)
    ## Theta function implementation with conditional
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2,cutoff, m, n)
        eps1,eps2 = get_dielectric(dielectric, k1, k2,cutoff, m, n)

        return 2.2*(((momentum^2 - k1^2 + k2^2)*(tanh(vel2*k2/(2.0*temp))))/(momentum^2*eps1) + ((momentum^2 + k1^2 - k2^2)*(tanh(vel1*k1/(2.0*temp))))/(momentum^2*eps2))/(2.0*pi*sqrt((k1+k2)^2 - momentum^2))
    end
end

function dielectric_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64,temp::Float64)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, cutoff, m, n)

        tanh1 = tanh(vel1*k1/(2.0*temp))
        tanh2 = tanh(vel2*k2/(2.0*temp))

        denomin = sqrt((k1+k2)^2 - momentum^2)

        if abs(k1*vel1 - k2*vel2)<1e-5 # the integrand has a 0/0 singular form. To regulate that we use the limiting value.
            return 4.4*(-denomin*tanh1+ 2*k1*k2*(tanh1 + tanh2)/denomin)/(pi*momentum*(k1*vel1 + k2*vel2))
        else
            return 4.4*(denomin*(k2*vel2*tanh1 - k1*vel1*tanh2)/(k1*vel1 - k2*vel2) + 2*k1*k2*(tanh1 + tanh2)/denomin)/(pi*momentum*(k1*vel1 + k2*vel2))
        end
    end
end

temps = [0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]

for temp in temps
    ## Boundary values initialisation
    velocity[:,m] .= 1.0
    dielectric[:,m] .= 1.0

    velocity_integrand(velocity,dielectric,momentum,cutoff,phi,m,n) = velocity_integrand_t(velocity,dielectric,momentum,cutoff,phi,m,n,temp)
    dielectric_integrand(velocity,dielectric,momentum,cutoff,phi,m,n) = dielectric_integrand_t(velocity,dielectric,momentum,cutoff,phi,m,n,temp)

    ## solving exact FRG using an user defined function from FRGn Package
    rg_procedure(velocity,dielectric,velocity_integrand, dielectric_integrand ,m,n)
end
