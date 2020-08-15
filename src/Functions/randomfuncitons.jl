@doc raw"""
contains some of the useful functions that I use frequently.
"""
module functions

using PyCall, Interpolations
# using JLD, HDF5

export Hilbert

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

function interp(xlist::StepRangeLen,array::AbstractArray)
    return scale(interpolate(array,BSpline(Quadratic(Reflect(OnCell())))),xlist)
end

end