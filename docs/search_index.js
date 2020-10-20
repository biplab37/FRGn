var documenterSearchIndex = {"docs":
[{"location":"submodules/samplefunctions/#Sample-Functions","page":"Sample Functions","title":"Sample Functions","text":"","category":"section"},{"location":"submodules/samplefunctions/","page":"Sample Functions","title":"Sample Functions","text":"Pages = [\"samplefunctions.md\"]","category":"page"},{"location":"submodules/samplefunctions/#Bauer","page":"Sample Functions","title":"Bauer","text":"","category":"section"},{"location":"submodules/samplefunctions/","page":"Sample Functions","title":"Sample Functions","text":"FRGn.Bauer.velocity_integrand\nFRGn.Bauer.dielectric_integrand","category":"page"},{"location":"submodules/samplefunctions/#FRGn.Bauer.velocity_integrand","page":"Sample Functions","title":"FRGn.Bauer.velocity_integrand","text":"velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)\n\nThis function returns the integrand of the FRG equation for the velocity renormlisation.\n\nArgs\n\ndielectric (Array) : Array containing the dielectric values\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\ni          (Int64) : index for the running cutoff\n\n\n\n\n\n","category":"function"},{"location":"submodules/samplefunctions/#FRGn.Bauer.dielectric_integrand","page":"Sample Functions","title":"FRGn.Bauer.dielectric_integrand","text":"dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)\n\nThis function returns the integrand of the FRG equation for the dielectric function renormlisation.\n\nArgs\n\nvelocity   (Array) : The velocity array\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\ni          (Int64) : index for the running cutoff\n\n\n\n\n\n","category":"function"},{"location":"submodules/samplefunctions/#BosonInShell","page":"Sample Functions","title":"BosonInShell","text":"","category":"section"},{"location":"submodules/samplefunctions/","page":"Sample Functions","title":"Sample Functions","text":"FRGn.BosonInShell.velocity_integrand\nFRGn.BosonInShell.dielectric_integrand","category":"page"},{"location":"submodules/samplefunctions/#FRGn.BosonInShell.velocity_integrand","page":"Sample Functions","title":"FRGn.BosonInShell.velocity_integrand","text":"velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)\n\nThis function returns the integrand of the FRG equation for the velocity renormlisation.\n\nArgs\n\nvelocity   (Array) : The velocity array\ndielectric (Array) : Array containing the dielectric values\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate (should run from 0 to pi/2)\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\n\n\n\n\n\n","category":"function"},{"location":"submodules/samplefunctions/#FRGn.BosonInShell.dielectric_integrand","page":"Sample Functions","title":"FRGn.BosonInShell.dielectric_integrand","text":"dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64)\n\nThis function returns the integrand of the FRG equation for the dielectric function renormlisation.\n\nArgs\n\nvelocity   (Array) : The velocity array\ndielectric (Array) : Array containing the dielectric values\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate (should run from 0 to pi/2)\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\n\n\n\n\n\n","category":"function"},{"location":"submodules/samplefunctions/#Finite-Temperature","page":"Sample Functions","title":"Finite Temperature","text":"","category":"section"},{"location":"submodules/samplefunctions/","page":"Sample Functions","title":"Sample Functions","text":"FRGn.FiniteTemp.velocity_integrand_t\nFRGn.FiniteTemp.dielectric_integrand_t","category":"page"},{"location":"submodules/samplefunctions/#FRGn.FiniteTemp.velocity_integrand_t","page":"Sample Functions","title":"FRGn.FiniteTemp.velocity_integrand_t","text":"velocity_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, temp::Float64)\n\nThis function returns the integrand of the FRG equation for the velocity renormlisation at a specific temperature.\n\nArgs\n\nvelocity   (Array) : The velocity array\ndielectric (Array) : Array containing the dielectric values\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate (should run from 0 to pi/2)\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\ni          (Int64) : index for the running cutoff\ntemp     (Float64) : temperature\n\n\n\n\n\n","category":"function"},{"location":"submodules/samplefunctions/#FRGn.FiniteTemp.dielectric_integrand_t","page":"Sample Functions","title":"FRGn.FiniteTemp.dielectric_integrand_t","text":"dielectric_integrand_t(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, temp::Float64)\n\nThis function returns the integrand of the FRG equation for the dielectric function renormlisation at a specific temperature.\n\nArgs\n\nvelocity   (Array) : The velocity array\ndielectric (Array) : Array containing the dielectric values\nmomentum (Float64) : momentum value\ncutoff   (FLoat64) : running cutoff\nphi      (Float64) : angular coordinate\nm          (Int64) : number of cutoffs\nn          (Int64) : number of momenta\ni          (Int64) : index for the running cutoff\ntemp     (Float64) : temperature\n\n\n\n\n\n","category":"function"},{"location":"usage/#User-Guide","page":"Usage","title":"User Guide","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Pages = [\"usage.md\"]","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"To be able  solve exact FRG equations one need to follow the steps below.","category":"page"},{"location":"usage/#Preparation","page":"Usage","title":"Preparation","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Before starting coding one need to do the necessary analytical calculations in order to obtain the exact FRG equations that one wants to solve. Then one need to simplify the equations such as to mould them into the following form","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"partial_Lambdaf_i = int_0^pi2hatf_i","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"where f_i's are the observables that we want to solve and hatf_i are the integrands corrosponding to the observables. Note the limit of the integration is from 0 to pi2. If one needs to solve a problem which can not be turned into this form please raise an issue. Once we have the closed set of equations we are ready to solve them.","category":"page"},{"location":"usage/#Initialisation","page":"Usage","title":"Initialisation","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Note that the observables f_i's are dependent on the momentum and the running cutoff. So we would like them to store in a 2-dimensional array. So the code would look something like","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"using FRGn\n\nm = 123 # no of cutoffs\nn = 456 # no of momenta\nobservable1 = zeros(n,m)\nobservable2 = zeros(n,m)\n...\n","category":"page"},{"location":"usage/#Define-integrands","page":"Usage","title":"Define integrands","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Once the observables are defined the integrands on the rhs of FRG equation have to be provided. This should be done as a function of the following form","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"function integrand1(observable1 ,observable2, momentum, cutoff, phi, m, n) \n\t...\nend","category":"page"},{"location":"usage/#Solving","page":"Usage","title":"Solving","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Before solving the differential equations we need to supply the boundary condition. Boundary condition are comes from the fact that we know the observables at the highest cutoff i.e. before starting the renormalisation group procedure. For example","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"observable1[:,m] .= 1","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"where we redefined the variable so that the value of observable was unity at the highest cutoff.","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"Finally, to solve the FRG equations we call the rg_procedure funciton","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"rg_procedure(observable1, observable2, integrand1, integrand2, m, n)","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"This will update the observables in place.","category":"page"},{"location":"usage/#Visualise-Results","page":"Usage","title":"Visualise Results","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"To visulaise the result one can plot the resultant observables for specific cases using Plots or PyPlot. Here we provide two simple function to plot the renormlised velocity and dielectric function as follows","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"plot_velocity(velocity[:,1],\"renormalised_velocity.png\")","category":"page"},{"location":"usage/","page":"Usage","title":"Usage","text":"this will also save the figure in the given name. Omiting that variable will not save the figure.","category":"page"},{"location":"usage/#Further-Examples","page":"Usage","title":"Further Examples","text":"","category":"section"},{"location":"usage/","page":"Usage","title":"Usage","text":"Further examples are given in the example section of the code.","category":"page"},{"location":"indices/#List-of-Functions","page":"Indices","title":"List of Functions","text":"","category":"section"},{"location":"indices/","page":"Indices","title":"Indices","text":"List of the functions available","category":"page"},{"location":"indices/","page":"Indices","title":"Indices","text":"Order = [:function]","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Further examples can be found in the example section of the code.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Bauer\nBoson In Shell\nFinite Temperature","category":"page"},{"location":"submodules/getveleps/#GetVelEps","page":"GetVelEps","title":"GetVelEps","text":"","category":"section"},{"location":"submodules/getveleps/","page":"GetVelEps","title":"GetVelEps","text":"get_velocity\nget_dielectric\nfetch_value","category":"page"},{"location":"submodules/getveleps/#FRGn.get_velocity","page":"GetVelEps","title":"FRGn.get_velocity","text":"get_velocity(velocity::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)\n\nReturns the velocities at momenta k1 and k2 at the previous cutoff (indexed m-i+2).\n\n\n\n\n\nget_velocity(velocity::Array{Float64,2}, k::Float64, cutoff::Float64, m::Int64, n::Int64)\n\nReturns the velocities at momenta k at the previous cutoff (indexed m-i+2).\n\n\n\n\n\n","category":"function"},{"location":"submodules/getveleps/#FRGn.get_dielectric","page":"GetVelEps","title":"FRGn.get_dielectric","text":"get_dielectric(dielctric::Array{Float64,2}, k1::Float64, k2::Float64, cutoff::Float64, m::Int64, n::Int64)\n\nReturns the dielectric functions at momenta k1 and k2 at the previous cutoff (indexed m-i+2).\n\n\n\n\n\nget_dielectric(dielctric::Array{Float64,2}, k::Float64,cutoff::Float64, m::Int64, n::Int64)\n\nReturns the dielectric functions at momenta k at the previous cutoff (indexed m-i+2).\n\n\n\n\n\n","category":"function"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"The install process is very easy","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Clone the repository","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"git clone https://github.com/biplab37/FRGn.git","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Add the file path to the PATH of julia.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"julia> push!(LOAD_PATH,location/to/the/file)","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"One should add this line to the file \"~/.julia/config/startup.jl\" to avoid running the command in each julia session.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"After that you can add this package","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"julia> ]add FRGn","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"To import or use the functions in the package","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"julia> using FRGn","category":"page"},{"location":"submodules/plotting/#Plotting","page":"Plotting","title":"Plotting","text":"","category":"section"},{"location":"submodules/plotting/","page":"Plotting","title":"Plotting","text":"Some predefined format of plotting to avoid labeling axes everytime you plot.","category":"page"},{"location":"submodules/plotting/","page":"Plotting","title":"Plotting","text":"plot_velocity\nplot_dielectric","category":"page"},{"location":"submodules/plotting/#FRGn.plot_velocity","page":"Plotting","title":"FRGn.plot_velocity","text":"plot_velocity(velocity[,name])\n\nPlots the velocity as a function of momentum. Automatically labels the graph. If name is given then saves the figure with the given name.\n\nExample\n\njulia> plot_velocity(velocity,\"renormalised_velocity.pdf\")\n\n\n\n\n\n","category":"function"},{"location":"submodules/plotting/#FRGn.plot_dielectric","page":"Plotting","title":"FRGn.plot_dielectric","text":"plot_dielectric(dielectric[,name])\n\nPlots the dielectric as a function of momentum and automatically labels the graph. If name is given then saves the figure with the given name.\n\nExample\n\njulia> plot_dielectric(dielectric,\"renormalised_dielectric.pdf\")\n\n\n\n\n\n","category":"function"},{"location":"submodules/rgprocedure/#RGProcedure","page":"RGProcedure","title":"RGProcedure","text":"","category":"section"},{"location":"submodules/rgprocedure/","page":"RGProcedure","title":"RGProcedure","text":"rg_procedure","category":"page"},{"location":"submodules/rgprocedure/#FRGn.rg_procedure","page":"RGProcedure","title":"FRGn.rg_procedure","text":"rg_procedure(velocity::Array{Float64,2}, dielectric::Array{Float64,2}, velocity_integrand::Function, dielectric_integrand::Function, m::Int64, n::Int64)\n\nThis function solves the exact FRG equations for velocity and dielectric function where one should input the functions velocity_integrand and dielectric_integrand which returns the integrands defined as below. Note that the angular coordinate phi should run from 0 to pi/2.  \n\n    velocity_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff:;Float64, phi:: Float64, m::Int64, n::Int64, i::Int64) \n    dielectric_integrand(velocity::Array{Float64,2},dielectric::Array{Float64,2}, momentum::Float64, cutoff:;Float64, phi:: Float64, m::Int64, n::Int64, i::Int64) \n\n\n\n\n\nrg_procedure(functions::Array{Array{Float64,2},1}, functions_integrand::Array{Function,1}, m::Int64, n::Int64)\n\nThe function solves the exact FRG equations for the given set of functions(for example velocity, dielectric) input as an array. Functions returning the integrand should be provided in a array too. The integrand functions should have the following from\n\n    function_integrand(functions::Array{Array{Float64,2},1}, momentum::Float64, cutoff::Float64, phi::Float64, m::Int64, n::Int64, i::Int64)\n\n\n\n\n\n","category":"function"},{"location":"#FRGn-Documentation","page":"Home","title":"FRGn Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A small package that contains some of the subroutines and function which I frequently use in my code to numerically solve the exact Funcitonal Rrenormalisation Group equation in Graphene. The equations are motivated from the paper by Bauer et al [1]. Here apart from reproducing the result by them we also extend it to finite temperature and introduce the frequency dependence in the interraction term. See the example page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[1]:  Carsten Bauer, Andreas Rückriegel, Anand Sharma, and Peter Kopietz. Non-perturbative renormalization group calculation of quasiparticle velocity and di-electric function of graphene. Phys. Rev. B, 92:121409, Sep 2015.","category":"page"},{"location":"#Indices","page":"Home","title":"Indices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"indices.md\"]","category":"page"}]
}