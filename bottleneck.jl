### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 162029e0-9c77-11eb-1937-69fb4af41087
begin
	using LinearAlgebra
	using LightGraphs
	using GraphRecipes
	using Plots
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 84fe8c10-3715-4c94-8a59-7aa87ba50ff7
pyplot()

# ╔═╡ d087a3f9-f88a-43ef-828e-e269039fdc45
function henrici(A::AbstractMatrix; normalize=true)
    eig = eigen(A)
    frob = norm(A) # by default, calling norm on a matrix returns the Frobenius norm
    if normalize
        return sqrt( max( frob^2 - sum(abs.(eig.values).^2), 0.0) )/frob
    else
        return sqrt( max( frob^2 - sum(abs.(eig.values).^2), 0.0) )
    end
end

# ╔═╡ 55e9a97f-e770-47ea-a28b-4aa671f163ff
function h_exact(n, ε; normalize = true)
	if normalize
		return sqrt(n-1 + ε^2 - n*ε^(2/n)) / sqrt(n-1+ε^2)
	else
		return sqrt(n-1 + ε^2 - n*ε^(2/n))
	end
end

# ╔═╡ 04c27d15-bc70-45de-93ac-eb4ef7c6de71
h_(n; normalize=true) = x -> h_exact(n,x, normalize=normalize)

# ╔═╡ 8fa1ff8d-2da3-4299-a072-b07e3323df65
@bind M Slider(5:100, show_value=true)

# ╔═╡ 782234ec-7ff1-4835-b3e3-d76aa32bfbc4
begin
	g = cycle_digraph(M)
	A = Array{Float64}(adjacency_matrix(g));
end

# ╔═╡ 3db7c5d9-f362-41af-b4f5-2dd8957ef981
begin
	local labels = Dict((5,1)=> L"\varepsilon")
	graphplot(g, 
		edgelabel=labels,
		fontsize=12,
		edgelabel_box=false,
		edgelabel_offset=0.05,
		method=:circular
	)
	#savefig("images/bottleneck-graph.pdf")
end

# ╔═╡ 122d15db-392b-42a9-ab5d-363d387728eb
begin
	local n = 50
	local ε = 10.0 .^ LinRange(-16, 0, n)
	local h = Vector(undef, n)
	for i in 1:n
		A[M,1] = ε[i]
		h[i] = henrici(A)
	end
	plot(ε, (h_(M)).(ε),
		label="Analytical",
		xlabel=L"\varepsilon",
		ylabel=L"\hat{h}",
		title="Normalised departure from normality for N =$(M)",
		xscale=:log10,
		#yscale=:log10
	)
	scatter!(ε, h, label="Base Formula")
end

# ╔═╡ 26bb659f-3c19-42b7-aa9e-b2d40adec17c
begin
	local n = 100
	local N = [10,50,100,500,1000]
	local ε = 10.0 .^ LinRange(-16, 0, n)
	local h = Array{Float64}(undef, length(N), n)
	for j in 1:n, i in 1:length(N)
		h[i,j] = h_exact(N[i], ε[j])
	end
	local p = plot(;
		xlabel=L"\varepsilon",
		ylabel=L"\hat{h}",
		xscale=:log10,
	)
	for i in 1:length(N)
		plot!(p, ε, h[i,:], label="N=$(N[i])")
	end
	#savefig(p, "images/bottleneck-henrici.pdf")
	p
end

# ╔═╡ 3e65e0d0-37fd-4a65-97f4-9709dad4d3c8
begin
	local ε = [-2.5, -5.0, -10.0, -15.0]
	local n = 100000
	local h = Array{Float64,2}(undef, length(ε), n)
	for j in 1:n, i in 1:length(ε)
		h[i,j] = h_exact(j+1, 10.0^ε[i])
	end
	local p = plot(;
		xlabel=L"N",
		ylabel=L"\hat{h}",
		yscale=:log10,
		xscale=:log10
	)
	for i in 1:length(ε)
		plot!(p, 2:(n+1), h[i,:], label=L"\varepsilon = 10^{%$(ε[i])}")
	end
	#savefig(p, "images/bottleneck-henrici-2.pdf")
	p
end

# ╔═╡ 645baff1-ee56-4aca-a367-f6e1ebedc661
begin
	local ε = 10.0 .^ -LinRange(1.0,16.0,100)
	local N = 10 .^ LinRange(1.0,16.0, 100)
	contour(N, ε, (x,y) -> log10(h_exact(x,y)),
		xscale=:log10,
		yscale=:log10,
		fill=true,
		xlabel=L"N",
		ylabel=L"\varepsilon",
		title=L"\log_{10} \hat{h}",
		levels=50,
		grid=false
	)
	#savefig("images/bottleneck-contour.pdf")
end

# ╔═╡ 00ec2f66-38b1-423a-a00c-3a91866e5582
begin
	local ε = 10.0 .^ -LinRange(1.0,16.0,100)
	local N = 10 .^ LinRange(1.0,16.0, 100)
	contour(N, ε, (x,y) -> log10(sqrt(x-1+y^2 - x*y^(2/x))),
		xscale=:log10,
		yscale=:log10,
		fill=true,
		xlabel=L"N",
		ylabel=L"\varepsilon",
		title=L"\log_{10} h",
		levels=50,
		grid=false
	)
end

# ╔═╡ 837ae74e-93f4-4430-b433-8be904b098ba
begin
	local ε = BigFloat.([-2.5, -5.0, -10.0, -15.0])
	local N = BigFloat.(10 .^ LinRange(1.0,20.0, 200))
	local h = Array{BigFloat,2}(undef, length(ε), length(N))
	f(x,y) = sqrt(x-1 + y^2 - x * y^(2/x))
	for j in 1:length(N), i in 1:length(ε)
		h[i,j] = f(N[j], 10.0^ε[i])
	end
	local p = plot(;
		xlabel=L"N",
		ylabel=L"h",
		#yscale=:log10,
		xscale=:log10
	)
	for i in 1:length(ε)
		plot!(p, N, h[i,:], label=L"\varepsilon = 10^{%$(ε[i])}")
	end
	#savefig(p, "images/bottleneck-henrici-3.pdf")
	p
end

# ╔═╡ Cell order:
# ╠═162029e0-9c77-11eb-1937-69fb4af41087
# ╠═84fe8c10-3715-4c94-8a59-7aa87ba50ff7
# ╠═d087a3f9-f88a-43ef-828e-e269039fdc45
# ╠═782234ec-7ff1-4835-b3e3-d76aa32bfbc4
# ╠═3db7c5d9-f362-41af-b4f5-2dd8957ef981
# ╠═55e9a97f-e770-47ea-a28b-4aa671f163ff
# ╠═04c27d15-bc70-45de-93ac-eb4ef7c6de71
# ╠═8fa1ff8d-2da3-4299-a072-b07e3323df65
# ╠═122d15db-392b-42a9-ab5d-363d387728eb
# ╠═26bb659f-3c19-42b7-aa9e-b2d40adec17c
# ╠═3e65e0d0-37fd-4a65-97f4-9709dad4d3c8
# ╠═645baff1-ee56-4aca-a367-f6e1ebedc661
# ╠═00ec2f66-38b1-423a-a00c-3a91866e5582
# ╠═837ae74e-93f4-4430-b433-8be904b098ba
