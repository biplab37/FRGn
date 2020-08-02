module GetVelEps

export get_velocity, get_dielectric, fetch_value

"""
    get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)

Returns the velocities at momenta k1 and k2 at the previous cutoff (indexed m-i+2).
"""
function get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64,n::Int64)

    vel1 = get_velocity(velocity,k1,cutoff,m,n)
    vel2 = get_velocity(velocity,k2,cutoff,m,n)

    return vel1,vel2
end

"""
    get_velocity(velocity::Array{Float64,2}, k::Float64, cutoff::Float64, m::Int64, n::Int64)

Returns the velocities at momenta k at the previous cutoff (indexed m-i+2).
"""
function get_velocity(velocity::Array{Float64,2}, k::Float64, cutoff::Float64, m::Int64, n::Int64)

    index = Int64(round(k*n))
    i = Int64(round(m + 1 - (cutoff*m)))

    i = (i==0 || i==1) ? 2 : i
    index = (index==0) ? 1 : index

    if index<n
        vel::Float64 = velocity[index,m-i+2]
    else
        vel = 1.0 # if the index goes out of the boundary take velocity to be v_F.
    end

    return vel
end

"""
    get_dielectric(dielctric::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)

Returns the dielectric functions at momenta k1 and k2 at the previous cutoff (indexed m-i+2).
"""
function get_dielectric(dielectric::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)

    eps1 = get_dielectric(dielectric,k1,cutoff,m,n)
    eps2 = get_dielectric(dielectric,k2,cutoff,m,n)

    return eps1,eps2
end

"""
    get_dielectric(dielctric::Array{Float64,2}, k::Float64,cutoff::Float64, m::Int64, n::Int64)

Returns the dielectric functions at momenta k at the previous cutoff (indexed m-i+2).
"""
function get_dielectric(dielectric::Array{Float64,2}, k::Float64, cutoff::Float64, m::Int64,n::Int64)

    index = Int64(round(k*n))
    i = Int64(round(m + 1 - (cutoff*m)))

    i = (i==0 || i==1) ? 2 : i
    index = (index==0) ? 1 : index

    if index<n
        eps1::Float64 = dielectric[index,m-i+2]
    else
        eps1 = 1.0 # if the index goes out of the boundary take dielectric to be the free space one.
    end

    return eps1
end

"""
    fetch_value(input, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)

This function returns the values of the input at given momenta k1 and k2 and the previous cutoff value.

"""

function fetch_value(input, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)

    value1 = fetch_value(input,k1,cutoff,m,n)
    value2 = fetch_value(input,k2,cutoff,m,n)

    return value1, value2
end

"""
    fetch_value(input, k::Float64, cutoff::Float64, m::Int64, n::Int64)

This function returns the value of the input at given momentum k and the previous cutoff value.

"""

function fetch_value(input, k::Float64, cutoff::Float64, m::Int64, n::Int64)
    index = Int64(round(k*n))
    i = Int64(round(m + 1 - (cutoff*m)))

    i = (i==0 || i==1) ? 2 : i
    index = (index==0) ? 1 : index

    if index<n
        value::Float64 = input[index,m-i+2]
    else
        value = 1.0 # if the index goes out of the boundary take dielectric to be the free space one.
    end

    return value
end

end