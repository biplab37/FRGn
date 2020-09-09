@doc raw"""
# RGProcedure
This module provides the main function which solves the exact FRG equations.
"""

using QuadGK

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

# Same function if you dont want to keep the cutoff dependence
function rg_procedure(velocity::Array{Float64,1}, dielectric::Array{Float64,1}, velocity_integrand::Function, dielectric_integrand::Function, m::Int64, n::Int64)

    dcutoff::Float64 = 1.0/m
    new_velocity = zeros(n)
    new_dielectric = zeros(n)

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi,m,n)
            dielectric_integrand_phi(phi) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi,m,n)
            ## Solving ODE's using Euler method
            new_velocity[j] = velocity[j] + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            new_dielectric[j] = dielectric[j] + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-4)[1]
        end
        velocity , dielectric = new_velocity, new_dielectric
    end
end

# for functional inputs
function rg_procedure(velocity::Function, dielectric::Function, velocity_integrand::Function, dielectric_integrand::Function, m::Int64,n::Int64 = 100)

    new_velocity = zeros(n)
    new_dielectric = zeros(n)
    dcutoff = 1/m
    for i in 2:m
        cutoff = Float64(m - i +1)/m
        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi)

            ## Solving ODE's using Euler method
            new_velocity[j] = velocity(momentum) + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            new_dielectric[j] = dielectric(momentum) + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity = interp(n,new_velocity)
        dielectric = interp(n,new_dielectric)
    end

    return velocity, dielectric
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

function rg_procedure(functions::Array{Array{Float64,1},1}, functions_integrand::Array{Function,1}, m::Int64, n::Int64)

    dcutoff::Float64 = 1.0/m

    if length(functions) != length(functions_integrand)
        error("Number of functions does not match number of integrands!")        
    end

    new_functions = functions

    for i in 2:m

        cutoff = Float64(m - i +1)/m

        for j in 1:n

            momentum = Float64(j)/n

            for k in 1:length(functions)
                integrand_phi(phi) = functions_integrand[k](functions,momentum,cutoff,phi,m,n)
                new_functions[k,j] = functions[k][j] + dcutoff*quadgk(integrand_phi,0.,pi/2.,rtol=1e-4)[1]
            end
        end
        functions = new_functions
    end
end

function rg_procedure(property::Array{Float64,2},integrand_function::Function,m::Int64,n::Int64)

    dcutoff::Float64 = 1.0/n

    for i in 2:m

        cutoff = Float64(m-i+1)/m

        for j in 1:n

            momentum = Float64(j)/n

            integrand_phi(phi) = integrand_function(property,momentum,cutoff,phi,m,n)
            property[j,m-i+1] = property[j,m-i+2] + dcutoff*quadgk(integrand_phi,0.0,pi/2,rtol=1e-4)[1]
            
        end
        
    end
    
end

function rg_procedure(property::Array{Float64,2},integrand_function::Function,m::Int64,n::Int64,extra_arguments::Array{Float64})

    dcutoff::Float64 = 1.0/n

    for i in 2:m

        cutoff = Float64(m-i+1)/m

        for j in 1:n

            momentum = Float64(j)/n

            integrand_phi(phi) = integrand_function(property,momentum,cutoff,phi,m,n,extra_arguments)
            property[j,m-i+1] = property[j,m-i+2] + dcutoff*quadgk(integrand_phi,0.0,pi/2,rtol=1e-4)[1]
            
        end
        
    end
    
end

function rg_procedure(property::Array{Float64,1},integrand_function::Function,momentum::Float64,m::Int64,n::Int64,extra_arguments::Array{Float64})

    dcutoff::Float64 = 1.0/n

    for i in 2:m

        cutoff = Float64(m-i+1)/m

        integrand_phi(phi) = function_integrand(property,momentum,cutoff,phi,m,n,extra_arguments)
        property[m-i+1] = property[m-i+2] + dcutoff*quadgk(integrand_phi,0.0,pi/2,rtol=1e-4)[1]
    
    end
    
end

function rg_procedure_rk2(velocity::Function, dielectric::Function, velocity_integrand::Function, dielectric_integrand::Function, m::Int64,n::Int64 = 100)

    new_velocity = zeros(n)
    new_dielectric = zeros(n)
    dcutoff = 1/m
    for i in 2:m
        cutoff = Float64(m - i +1)/m
        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi)

            ## Solving ODE's using Euler method
            new_velocity[j] = velocity(momentum) + 1/2*dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            new_dielectric[j] = dielectric(momentum) + 1/2*dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        mid_velocity = interp(n,new_velocity)
        mid_dielectric = interp(n,new_dielectric)

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(mid_velocity,mid_dielectric,momentum,cutoff + 1/2*dcutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(mid_velocity,mid_dielectric,momentum,cutoff + 1/2*dcutoff,phi)

            ## Solving ODE's using Euler method
            new_velocity[j] = velocity(momentum) + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            new_dielectric[j] = dielectric(momentum) + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity = interp(n,new_velocity)
        dielectric = interp(n,new_dielectric)
    end

    return velocity, dielectric
end

function rg_procedure_rk4(velocity::Function, dielectric::Function, velocity_integrand::Function, dielectric_integrand::Function, m::Int64,n::Int64 = 100)

    temp_velocity = zeros(n)
    temp_dielectric = zeros(n)
    dcutoff = 1/m
    for i in 2:m
        cutoff = Float64(m - i +1)/m
        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity,dielectric,momentum,cutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity,dielectric,momentum,cutoff,phi)

            ## Solving ODE's using Euler method
            temp_velocity[j] = velocity(momentum) + 1/2*dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            temp_dielectric[j] = dielectric(momentum) +1/2*dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity_1 = interp(n,temp_velocity)
        dielectric_1 = interp(n,temp_dielectric)

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity_1,dielectric_1,momentum,cutoff + 1/2*dcutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity_1,dielectric_1,momentum,cutoff + 1/2*dcutoff,phi)

            ## Solving ODE's using Euler method
            temp_velocity[j] = velocity(momentum) + 1/2*dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            temp_dielectric[j] = dielectric(momentum) +1/2*dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity_2 = interp(n,temp_velocity)
        dielectric_2 = interp(n,temp_dielectric)

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity_2,dielectric_2,momentum,cutoff + 1/2*dcutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity_2,dielectric_2,momentum,cutoff + 1/2*dcutoff,phi)

            ## Solving ODE's using Euler method
            temp_velocity[j] = velocity(momentum) + dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            temp_dielectric[j] = dielectric(momentum) + dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity_3 = interp(n,temp_velocity)
        dielectric_3 = interp(n,temp_dielectric)

        for j in 1:n

            momentum = Float64(j)/n

            velocity_integrand_phi(phi::Float64) = velocity_integrand(velocity_3,dielectric_3,momentum,cutoff+dcutoff,phi)
            dielectric_integrand_phi(phi::Float64) = dielectric_integrand(velocity_3,dielectric_3,momentum,cutoff+dcutoff,phi)

            ## Solving ODE's using Euler method
            temp_velocity[j] =  velocity(momentum) - 1/2*dcutoff*quadgk(velocity_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
            temp_dielectric[j] =  dielectric(momentum) - 1/2*dcutoff*quadgk(dielectric_integrand_phi,0.,pi/2.,rtol=1e-3)[1]
        end
        velocity_4 = interp(n,temp_velocity)
        dielectric_4 = interp(n,temp_dielectric)

        velocity_5(x) = 1/3*(velocity_1(x) + 2*velocity_2(x) + velocity_3(x) - velocity_4(x))
        dielectric_5(x) = 1/3*(dielectric_1(x) + 2*dielectric_2(x) + dielectric_3(x) - dielectric_4(x))
        velocity,dielectric = velocity_5,dielectric_5
    end

    return velocity, dielectric
end