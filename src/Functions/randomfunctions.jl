@doc raw"""
contains some of the useful functions that I use frequently.
"""
module functions

using PyCall
# using JLD, HDF5

export Hilbert, interp

@doc raw"""
    Hilbert(list::Array{Float,1})

This function return the hilbert transformation of a given list
"""
function Hilbert(list::Array{Float64,1})

    ft = pyimport("scipy.fftpack")
    
    return ft.hilbert(list)
end

@doc raw"""
    get_colormap(cmap::String)

This function returns the colormap from matplotlib given a string for the name of the colormap. Default is "coolwarm"
"""
function get_colormap(cmap::String = "coolwarm")

	cm = pyimport("matplotlib.cm")

	return cm.get_cmap(cmap)
end

# function load(filename::String, name::String)

# 	extension = split(filename,".")[2]

# 	if extension == "jld"
# 		return load(filename, name)
# 	elseif extension =="h5"
# 		return h5read(filename, name)
# 	end
# end

function itp(momentum::Float64,n::Int64,list::Array)

    index1 = Int64(floor(n*momentum))
    index2 = Int64(ceil(n*momentum))

    if index1==0
        return list[1]
    elseif index2>n
        return list[n]
    else
        if index1 == index2
            return list[index1]
        else
            return list[index1] + (list[index2] - list[index1])*(n*momentum - index1)
        end
    end
end

@doc raw"""
    function interp(n::Int64,list::Array)
This function returns the interpolated function given an Array.
"""
function interp(n::Int64,list::Array)
    itp2(momentum::Float64) = itp(momentum,n,list)
    return itp2
end

end