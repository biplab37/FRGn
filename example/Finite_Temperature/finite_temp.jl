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

using PyPlot
using PyCall

matplotlib = pyimport("matplotlib")
cmap = matplotlib.cm.get_cmap("coolwarm") #matplotlib colormap

fig1 = figure()
fig2 = figure()

ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

# using HDF5
temps = [0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]

for temp in temps
    ## Boundary values initialisation
    velocity[:,m] .= 1.0
    dielectric[:,m] .= 1.0

    velocity_integrand(velocity,dielectric,momentum,cutoff,phi,m,n) = velocity_integrand_t(velocity,dielectric,momentum,cutoff,phi,m,n,temp)
    dielectric_integrand(velocity,dielectric,momentum,cutoff,phi,m,n) = dielectric_integrand_t(velocity,dielectric,momentum,cutoff,phi,m,n,temp)

    ## solving exact FRG using an user defined function from FRGn Package
    rg_procedure(velocity,dielectric,velocity_integrand, dielectric_integrand ,m,n)



    ## Plots using user defined functions in FRGn Package

    ax1.plot(range(0,stop=1,length=n),velocity[:,1],label="t = $temp",color=cmap(10*temp))

    ax2.plot(range(0,stop=1,length=n),dielectric[:,1],label="t = $temp",color=cmap(10*temp))

    # h5write("finite.h5",string("t",temp,"/velocity"),velocity)
    # h5write("finite.h5",string("t",temp,"/dielectric"),dielectric)
end

ax1.set_title("Temperature Dependence of Renormlised Velocity")
ax1.set_xlabel(L"$k/\Lambda_0$")
ax1.set_ylabel(L"$\dfrac{v_{\Lambda \to 0}(k)}{v_F}$")
ax1.set_xlim(0,1)
ax1.legend()
fig1.savefig("temperature_dependence_velocity.pdf")

ax2.set_title("Temperature Dependence of Renormalised Dielectric Function")
ax2.set_xlabel(L"$q/\Lambda_0$")
ax2.set_ylabel(L"$\epsilon_{\Lambda \to 0}(q)$")
ax2.set_ylim(1,2.45)
ax2.set_xlim(0,1)
ax2.legend()
fig2.savefig("temperature_dependence_dielectric.pdf")
