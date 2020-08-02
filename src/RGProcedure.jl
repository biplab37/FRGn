@doc raw"""
# RGProcedure
This module provides the main function which solves the exact FRG equations.
"""
module RGProcedure

using QuadGK

export rg_procedure

@doc raw"""
    rg_procedure(velocity::Array{Float64,2}, dielectric::Array{Float64,2}, velocity_integrand::Function, dielectric_integrand::Function, m::Int64, n::Int64)
This function solves the exact FRG equations for velocity and dielectric function where one should input the functions *velocity_integrand* and *dielectric_integrand* which returns the integrands defined as below. Note that the angular coordinate *phi* should run from 0 to pi/2.  
```julia
    velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff:;Float64, phi:: Float64, m::Int64, n::Int64, i::Int64) 
    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff:;Float64, phi:: Float64, m::Int64, n::Int64, i::Int64) 
```
"""
function rg_procedure(velocity::Array{Float64,2}, dielectric::Array{Float64,2}, velocity_integrand::Function, dielectric_integrand::Function, m::Int64, n::Int64)

	dcutoff::Float64 = 1.0/m

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi,m,n)
            dielectric_integrand_phi(phi) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi,m,n)
            ## Solving ODE's using Euler method
            velocity[j,m-i+1] = velocity[j,m-i+2] + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            dielectric[j,m-i+1] = dielectric[j,m-i+2] + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
        end
    end
end

@doc raw"""
    rg_procedure(functions::Array{Array{Float64,2},1}, functions_integrand::Array{Function,1}, m::Int64, n::Int64)
The function solves the exact FRG equations for the given set of functions(for example velocity, dielectric) input as an array. Functions returning the integrand should be provided in a array too. The integrand functions should have the following from
```julia
    function_integrand(functions::Array{Array{Float64,2},1}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)
```

"""
function rg_procedure(functions::Array{Array{Float64,2},1}, functions_integrand::Array{Function,1}, m::Int64, n::Int64)

    dcutoff::Float64 = 1.0/m

    if length(functions) != length(functions_integrand)
        error("Number of functions does not match number of integrands!")        
    end

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            for k in 1:length(functions)
                integrand_phi(phi) = functions_integrand[k](functions,momentum,cutoff,phi,m,n)
                functions[k][j,m-i+1] = functions[k][j,m-i+2] + dcutoff*quadgk(integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            end
        end
    end
end
end