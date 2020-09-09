### A Pluto.jl notebook ###
# v0.11.10

using Markdown
using InteractiveUtils

# ╔═╡ f0259856-eea2-11ea-128a-07dd90f42f90
using JLD, PyPlot, PyCall

# ╔═╡ be657e90-eea1-11ea-3d58-d9f43ecfc833
md" # Dielectric"

# ╔═╡ 0cf5508a-eea2-11ea-2482-97a055eb7e96
md"In the small frequency limit we can write the dielectric function as:
 
$\epsilon(i\omega,q) = X(q) + \omega^2 Y(q)$

where $X(q)$ and $Y(q)$ can be calculated numerically (see the [code](correction.jl)).
"

# ╔═╡ f8e03848-eea2-11ea-1ccf-4955fd01195c
dielectric,dielectric2= load("correction.jld","dielectric","dielectric2");

# ╔═╡ 37fab3a0-eea3-11ea-1c56-2b4e60c00d55
begin
	X = dielectric[:,1];
	Y = dielectric2[:,1];
end

# ╔═╡ 4e925b18-eea3-11ea-2f6c-fb66a7eccb98
md"*Pole correction approximation:*

$\frac{1}{\epsilon(\omega,q)} = 1 + \frac{A_q}{\omega^2 - \omega_q^2}$

Equating this with the previous equations in the small frequency limit we get

$A_q = - \frac{(X(q)-1)^2}{Y(q)} \quad \omega_q^2 = - \frac{X(q)(X(q) - 1)}{Y(q)}$

"

# ╔═╡ d308acd0-eea3-11ea-3986-a192a415edc6
begin
	Aq = - (X .- 1.0).^2 ./Y;
	ωq = -(X .* (X .- 1.0)) ./Y;
end

# ╔═╡ fa6fcd94-eea3-11ea-1175-7f94144235cb
begin
	close("all")
	plot(range(0,stop=1,length=length(Aq)),Aq)
	title(L"Renormalised $A_q$")
	xlabel(L"$\dfrac{k}{Λ_0}$")
	ylabel(L"$A_q$")
	gcf()
end

# ╔═╡ 6dd6c058-eea4-11ea-0b87-87b66cd6f395
begin
	close("all")
	plot(range(0,stop=1,length=length(ωq)),ωq)
	title(L"Renormalised $ω_q$")
	xlabel(L"$\dfrac{k}{Λ_0}$")
	ylabel(L"$ω_q$")
	gcf()
end

# ╔═╡ 962946b6-eea4-11ea-3787-857a5e4944ba
begin
	m = length(dielectric[1,:])
	n = length(dielectric[:,1])
	mpl = pyimport("matplotlib")
	cmap = mpl.cm.get_cmap("coolwarm")
	
	fig = figure(figsize=[10,12])
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)
	
	for i in 1:34:m
		Xp = dielectric[:,i]
		Yp = dielectric2[:,i]
		Aqp = -(Xp .- 1.0).^2 ./Yp;
		ωqp = -(Xp .* (Xp .- 1.0)) ./Yp;
		ax1.plot(range(0,stop=1,length=n),Aqp,color=cmap((i)/m),label="Λ = $(round((i)/m,digits=1)) Λ₀")
		ax2.plot(range(0,stop=1,length=n),ωqp,color=cmap((i)/m),label="Λ = $(round((i)/m,digits=1)) Λ₀")
	end
	ax1.legend()
	ax1.set_title(L"Change in $A_q$ in the renormalisation process")
	ax1.set_xlabel(L"$\dfrac{k}{\Lambda_0}$")
	ax1.set_ylabel(L"$A_{\Lambda}(q)$")
	ax2.set_title(L"Change in $\omega_q$ in the renormalisation process")
	ax2.set_xlabel(L"$\dfrac{k}{\Lambda_0}$")
	ax2.set_ylabel(L"$\omega_{\Lambda}(q)$")
	ax2.legend();

	gcf()
end

# ╔═╡ Cell order:
# ╟─be657e90-eea1-11ea-3d58-d9f43ecfc833
# ╟─0cf5508a-eea2-11ea-2482-97a055eb7e96
# ╠═f0259856-eea2-11ea-128a-07dd90f42f90
# ╠═f8e03848-eea2-11ea-1ccf-4955fd01195c
# ╠═37fab3a0-eea3-11ea-1c56-2b4e60c00d55
# ╟─4e925b18-eea3-11ea-2f6c-fb66a7eccb98
# ╠═d308acd0-eea3-11ea-3986-a192a415edc6
# ╟─fa6fcd94-eea3-11ea-1175-7f94144235cb
# ╟─6dd6c058-eea4-11ea-0b87-87b66cd6f395
# ╟─962946b6-eea4-11ea-3787-857a5e4944ba
