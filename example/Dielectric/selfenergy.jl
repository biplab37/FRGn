### A Pluto.jl notebook ###
# v0.11.10

using Markdown
using InteractiveUtils

# ╔═╡ 2503d004-eea7-11ea-3f0e-45c0924b425b
using JLD, PyPlot

# ╔═╡ 89bcd346-eea6-11ea-0e94-ddfa7040d71e
md"# Self Energy

###### Real and Imaginary part of self enrgy at k=0"

# ╔═╡ 32954eda-eea7-11ea-107f-b7f78a2afce7
velocity, Aq, ωq = load("renormalised_data.jld","velocity","Aq","ωq");

# ╔═╡ 44471034-eea7-11ea-2785-f13875bcd69d
n = length(velocity);

# ╔═╡ be1753cc-eea7-11ea-085a-a1e0951162a5
list = Float64[];

# ╔═╡ 4b0bc746-eea7-11ea-14e7-979c978efaef
for j in 1:10000
	omega = 10*j/10000 - 5
	int = 0.0
	for i in 1:n
		int += Aq[i]*((omega - ωq[i])/((i/n)^2*velocity[i]^2 + (omega - ωq[i])^2) + (omega + ωq[i])/((i/n)^2*velocity[i]^2 + (omega + ωq[i])^2))/(8*pi*ωq[i]^2)
	end
	int = int/n
	append!(list,int)
end

# ╔═╡ 0ee2fe1e-eea8-11ea-0569-a3b1836664d9
begin
	close("all")
	plot(range(-5,stop=5,length = 10000),list)
	title("Real part of the Self Energy at k=0");
	xlabel(L"\omega");
	ylabel(L"Re[\Sigma_{pa}(ω,0)]");
	gcf()
end

# ╔═╡ 68fc9a88-eea8-11ea-22be-b3f70d20a545
begin
	list3 = []
	ϵ = 1e-3
	close("all")
	for j in 1:10000
		omega = 10*j/10000 - 5
		int = 0.0
		for i in 1:n
		int += pi*omega*Aq[i]*(exp(-(omega + (i/n)*velocity[i]-ωq[i])^2/(4*ϵ)) - exp(-(omega - (i/n)*velocity[i] + ωq[i])^2/(4*ϵ)))/(8*(i/n)*sqrt(pi*ϵ)*ωq[i]*(ωq[i] - (i/n)*velocity[i]))
		end
		int = int/n
		append!(list3,int)
	end
	plot(range(-5,stop=5,length = 10000),list3)
	title("Imaginary part of the Self Energy at k=0");
	xlabel(L"\omega");
	ylabel(L"Im[\Sigma_{pa}(\omega,0)]");
	gcf()
end

# ╔═╡ Cell order:
# ╟─89bcd346-eea6-11ea-0e94-ddfa7040d71e
# ╠═2503d004-eea7-11ea-3f0e-45c0924b425b
# ╠═32954eda-eea7-11ea-107f-b7f78a2afce7
# ╠═44471034-eea7-11ea-2785-f13875bcd69d
# ╠═be1753cc-eea7-11ea-085a-a1e0951162a5
# ╠═4b0bc746-eea7-11ea-14e7-979c978efaef
# ╠═0ee2fe1e-eea8-11ea-0569-a3b1836664d9
# ╠═68fc9a88-eea8-11ea-22be-b3f70d20a545
