@doc raw"""
# FRGn package

This package contain the subroutines that are frequently shared between different codes to solve exact FRG equations. 

Following functions are available in this package.
 - get_velocity
 - get_dielectric
 - rg_procedure
 - plot_velocity
 - plot_dielctric

To know more about each of the functions use help. 
Typing '?' in the julia REPL will get you to help and then type the function name that you want help with.

```julia-repl
julia> ?
```

```julia-repl
help?> get_velocity
```
search: get_velocity

    get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, m::Int64, n::Int64, i::Int64)

Returns the velocities at momenta k1 and k2 at the previous cutoff (indexed m-i+2).
"""
module FRGn

## includes all the submodules
include("Plotting.jl")
include("GetVelEps.jl")
include("RGProcedure.jl")

## import all the submodules
using .Plotting
using .GetVelEps
using .RGProcedure

## make the functions global
export plot_velocity, plot_dielectric
export get_velocity, get_dielectric, get_velocityo, get_dielectrico
export rg_procedure

end # module
