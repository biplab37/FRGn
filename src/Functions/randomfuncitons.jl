module functions

using PyCall

export Hilbert
@doc raw"""
    Hilbert(list::Array{Float,1})

This function return the hilbert transformation of a given list
"""
function Hilbert(list::Array{Float64,1})
    ft = pyimport("scipy.fftpack")
    return ft.hilbert(list)
end

end