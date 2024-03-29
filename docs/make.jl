push!(LOAD_PATH,"../src/")
using Documenter, FRGn


makedocs(modules=[FRGn],
	sitename="FRGn",
	authors="Biplab Mahato",
	pages = Any[
		"Home"=>"index.md",
		"Installation"=>"installation.md",
		"Usage"=>"usage.md",
		"Examples"=>"examples.md",
		#hide("SubModules"=>["submodules/getveleps.md","submodules/plotting.md","submodules/rgprocedure.md"]),
		hide("Indices"=>"indices.md"),
		hide("Bauer"=>"bauer.md")])

deploydocs(repo="github.com/biplab37/FRGn.git")
