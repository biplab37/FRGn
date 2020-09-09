### A Pluto.jl notebook ###
# v0.11.10

using Markdown
using InteractiveUtils

# ╔═╡ 4e1d276e-ef7a-11ea-0a72-77f6bbfe3c11
push!(LOAD_PATH,"/media/biplab/e271932f-bbce-422e-967e-20895f7111c9/Code/julia/thesis");

# ╔═╡ 085dbd42-ef7a-11ea-275b-07a061e84ae0
using FRGn

# ╔═╡ 95f0997c-ef7d-11ea-1a66-4b3541557ea6
using PyPlot

# ╔═╡ 91d504dc-ef7f-11ea-2c84-b7390ddb8a17
md"# Bauer"

# ╔═╡ 4e9c530c-ef7f-11ea-01dd-5b7528b92c37
md"Load the path to the package"

# ╔═╡ 57ad6256-ef7f-11ea-1315-5bc573188872
md"import the package"

# ╔═╡ 5f657df8-ef7f-11ea-22c7-7b2eafa40617
md"Define constants"

# ╔═╡ 680d2c10-ef7c-11ea-0b9b-a72512c29cf7
begin
	const m = 313
	const n = 400;
end

# ╔═╡ 6938cf58-ef7f-11ea-1af2-b58f462150fe
md"Initialise the arrays"

# ╔═╡ 82d4940c-ef7c-11ea-04be-47c1566a192d
begin
	velocity = zeros(n,m);
	dielectric= zeros(n,m);
end

# ╔═╡ 7761833e-ef7f-11ea-1419-dfc0b63a67a7
md"Set up Boundary condition"

# ╔═╡ 51951b42-ef7d-11ea-13a3-e9b34bc84ebf
begin
	velocity[:,m] .= 1.0
	dielectric[:,m] .= 1.0
end

# ╔═╡ 5ea7ca76-ef7d-11ea-3dbb-cff21cadb6df
rg_procedure(velocity,dielectric,FRGn.Bauer.velocity_integrand, FRGn.Bauer.dielectric_integrand ,m,n)

# ╔═╡ 86b38c88-ef7f-11ea-1881-a9c305af5f26
md"###### Results"

# ╔═╡ 6fd36b20-ef7d-11ea-1945-61a47cda7d7a
begin
	close("all")
	plot(velocity[:,1])
	title("Renormalised Velocity")
	xlabel("k/Λ₀")
	ylabel(L"\dfrac{v_{\Lambda \to 0}(k)}{v_F}")
	gcf()
end

# ╔═╡ 0748dd78-ef7e-11ea-36be-e7ec10a44ad1
begin
	close("all")
	plot(dielectric[:,1])
	title("Renormalised Dielectric Function")
	xlabel("q/Λ₀")
	ylabel(L"ϵ_{Λ \to 0}(q)")
	gcf()
end

# ╔═╡ 3f186146-ef7f-11ea-16ba-0d5012f5a5fe
md"## Functional Form"

# ╔═╡ 53ed51e0-ef7e-11ea-2aac-4d3b7244c5a4
md"The results get significantly better if one uses interpolation functions at each steps of the RG process. This can be done by iniitialising the velocity and dielectric as funcitons"

# ╔═╡ 970d5ada-ef7e-11ea-0731-59e58577107e
begin
	initial_velocity(momentum) = 1.0
	initial_dielectric(momentum) = 1.0
end

# ╔═╡ acec51ce-ef7e-11ea-29cf-25e2a216f771
final_velocity, final_dielectric = rg_procedure(initial_velocity,initial_dielectric,FRGn.Bauer.velocity_integrand,FRGn.Bauer.dielectric_integrand,m,n);

# ╔═╡ fdc3205c-ef7e-11ea-3519-b713f219d5ca
begin
	close("all")
	plot(0:0.01:1,final_velocity.(0:0.01:1))
	title("Renormalised Velocity")
	xlabel("k/Λ₀")
	ylabel(L"\dfrac{v_{\Lambda \to 0}(k)}{v_F}")
	gcf()
end

# ╔═╡ 213b7f3c-ef7f-11ea-1c82-57e4e8532adc
begin
	close("all")
	plot(0:0.01:1,final_dielectric.(0:0.01:1))
	title("Renormalised Dielectric Function")
	xlabel("q/Λ₀")
	ylabel(L"ϵ_{Λ \to 0}(q)")
	gcf()
end

# ╔═╡ Cell order:
# ╟─91d504dc-ef7f-11ea-2c84-b7390ddb8a17
# ╟─4e9c530c-ef7f-11ea-01dd-5b7528b92c37
# ╠═4e1d276e-ef7a-11ea-0a72-77f6bbfe3c11
# ╟─57ad6256-ef7f-11ea-1315-5bc573188872
# ╠═085dbd42-ef7a-11ea-275b-07a061e84ae0
# ╟─5f657df8-ef7f-11ea-22c7-7b2eafa40617
# ╠═680d2c10-ef7c-11ea-0b9b-a72512c29cf7
# ╟─6938cf58-ef7f-11ea-1af2-b58f462150fe
# ╠═82d4940c-ef7c-11ea-04be-47c1566a192d
# ╟─7761833e-ef7f-11ea-1419-dfc0b63a67a7
# ╠═51951b42-ef7d-11ea-13a3-e9b34bc84ebf
# ╠═5ea7ca76-ef7d-11ea-3dbb-cff21cadb6df
# ╟─86b38c88-ef7f-11ea-1881-a9c305af5f26
# ╠═95f0997c-ef7d-11ea-1a66-4b3541557ea6
# ╠═6fd36b20-ef7d-11ea-1945-61a47cda7d7a
# ╠═0748dd78-ef7e-11ea-36be-e7ec10a44ad1
# ╟─3f186146-ef7f-11ea-16ba-0d5012f5a5fe
# ╟─53ed51e0-ef7e-11ea-2aac-4d3b7244c5a4
# ╠═970d5ada-ef7e-11ea-0731-59e58577107e
# ╠═acec51ce-ef7e-11ea-29cf-25e2a216f771
# ╠═fdc3205c-ef7e-11ea-3519-b713f219d5ca
# ╠═213b7f3c-ef7f-11ea-1c82-57e4e8532adc
