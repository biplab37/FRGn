module FiniteTemp

export velocity_integrand_t, dielectric_integrand_t

@doc raw"""
    velocity_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, temp::Float64)

This function returns the integrand of the FRG equation for the velocity renormlisation at a specific temperature.
## Args
    velocity   (Array) : The velocity array
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate (should run from 0 to pi/2)
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
    temp     (Float64) : temperature
"""
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

@doc raw"""
    dielectric_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, temp::Float64)

This function returns the integrand of the FRG equation for the dielectric function renormlisation at a specific temperature.
## Args
    velocity   (Array) : The velocity array
    dielectric (Array) : Array containing the dielectric values
    momentum (Float64) : momentum value
    cutoff   (FLoat64) : running cutoff
    phi      (Float64) : angular coordinate
    m          (Int64) : number of cutoffs
    n          (Int64) : number of momenta
    i          (Int64) : index for the running cutoff
    temp     (Float64) : temperature
"""
function dielectric_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64,temp::Float64)
    ## Theta function implementation
    if cos(phi)<=1 - 2*cutoff/momentum
        return 0.0
    else
        k1 = cutoff
        k2 = cutoff + cos(phi)*momentum

        vel1, vel2 = get_velocity(velocity, k1, k2, cutoff, m, n)
        # eps1,eps2 = get_dielectric(dielectric, k1, k2, m, n, i)

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

end