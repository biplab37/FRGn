# Installation

Follow the steps below to install the package

 - Clone the repository
```bash
git clone https://github.com/biplab37/FRGn.git
```

 - Add the file path to the PATH of julia.
```julia-repl
julia> push!(LOAD_PATH,location/to/the/file)
```
One should add this line to the file `~/.julia/config/startup.jl` to avoid running the command in each Julia session.

 - After that you can add this package
```julia-repl
julia> ]add FRGn
```
 - To import or use the functions in the package
```julia-repl
julia> using FRGn
```
