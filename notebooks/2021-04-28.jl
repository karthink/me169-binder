### A Pluto.jl notebook ###
# v0.14.3

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

# ╔═╡ d21e963c-a86a-11eb-3082-41cb8d8b1c17
begin
	using Plots, Pluto, PlutoUI, LinearAlgebra, Printf
end

# ╔═╡ c44dff8f-8e7a-4904-9a74-2cbf39250a03
md"""
## Plot of the vector field of $\dot{x} = A x$
"""

# ╔═╡ 8dc480f7-d8be-4826-bc0b-187ec47cf14d
md"Enter your $A$ matrix below:"

# ╔═╡ 2a237878-db73-4f29-87ea-cbcc49b642d0
md"""
 ``a_{1,1} = `` $(@bind a11 html"<input type=textbox size=4 placeholder='1'>") ``a_{1,2} = `` $(@bind a12 html"<input type=textbox size=4 placeholder='0'>")

 ``a_{2,1} = `` $(@bind a21 html"<input type=textbox size=4 placeholder='0'>") ``a_{2,2} = `` $(@bind a22 html"<input type=textbox size=4 placeholder='1'>")
"""

# ╔═╡ e75765e3-85fa-4eb1-abb0-7d740ec091cb
begin
	A = tryparse.(Float64,[ a11 a12; a21 a22 ])
end

# ╔═╡ 1c12a305-572d-44e5-8f25-7b0ecfc1313b
md"""### Vector field for $A$"""

# ╔═╡ bfac1751-6987-431b-967e-3305df530308
md"""
Normalize arrows? $(@bind n html"<input type=checkbox>")
"""

# ╔═╡ 86f58406-08c7-4adb-95ae-524c0d4666e7
md"arrow scale = $(@bind scale Slider(0.001:0.005:0.5, show_value=true, default=0.2))"

# ╔═╡ 08a1e247-4ba9-4693-bf11-06cc519320fc
md"""
## Jordan blocks
"""

# ╔═╡ 69451070-2f2a-4052-8b97-40625819db00
md"The matrix $J$ is a Jordan block (only one eigenvector) for $\epsilon = 0$. We investigate what happens to the eigenvectors and vector field as $\epsilon \to 0$."

# ╔═╡ 5479c733-f7c7-4ec7-8b4c-baff5fba9561
md"""
```math
J = \begin{bmatrix} -1 & 1 \\ \epsilon & -1 \end{bmatrix}
```
"""

# ╔═╡ 67adecf7-21e2-43aa-9dea-f881259ed939
md"arrow scale = $(@bind scaleJ Slider(0.001:0.005:0.5, show_value=true, default=0.2))"

# ╔═╡ 2b4e0633-025a-4ec9-8a01-24e11e1e8fd5
md"""
Normalize arrows? $(@bind nJ html"<input type=checkbox>")
"""

# ╔═╡ 1d59a29e-a930-481f-8093-3a232a1b6962
md"Change ϵ: $(@bind pow Slider(-14:1:4, show_value=true, default=0))"

# ╔═╡ 54d73ea3-86fe-4f7b-958a-6090e571bea2
begin
	ϵ = 2.0^pow
	md"#### ϵ = $(ϵ)"
end

# ╔═╡ f5a77a25-cc64-4d08-93fd-42beb3007cb4
md"Ignore"

# ╔═╡ fe63c64e-fef5-4abc-9636-d1413d40fc1b
md"### Helper functions"

# ╔═╡ a3e0aee0-b09c-4980-85b8-066cf9881757
begin
	J(ϵ) = [ -1 4; ϵ -1 ]
end

# ╔═╡ 8f29d603-c929-4004-a1f9-afd59e66b72f
function plotvec(plt, vecs)
	vecs[:,1] *= sign(vecs[1,1])
	vecs[:,2] *= sign(vecs[1,2])
	plot!(plt, t -> t * vecs[1,1], t -> t * vecs[2,1], 0, 1, lc=:teal, lw=2.5, arrow=true, label=false)
	plot!(plt, t -> t * vecs[1,2], t -> t * vecs[2,2], 0, 1, lc=:red, lw=2.5, arrow=true, label=false)
end

# ╔═╡ c9edac8f-8b56-4402-9315-9421f2bfbc7f
vectorfield(A) = vectorfield(A, identity, 0.2)

# ╔═╡ 8edbf3a3-6f75-4eff-887f-07005b860f43
function vectorfield(A, norm, scale)
	plt = plot(framestyle=:origin, aspect_ratio=1.0, xlims = (-1.5,1.5), ylims = (-0.75,0.75))
	plotvec(plt,eigvecs(A))
	gridnum = 25

	x = range(-2, stop=2, length=gridnum)
	y = range(-1, stop=1, length=gridnum)
	f(x, y) = norm(A * [x;y]) * scale
	quiver!(repeat(x,gridnum), vec(repeat(y',gridnum)), quiver=f, c=:blue, arrowscale=0.1, alpha=0.5, lc="#1d1f21")
	
	l1, l2 = eigvals(A)
	plot!(plt, title="λ = $(@sprintf "%.3f" l1), $(@sprintf "%.3f" l2)")
	return plt
end

# ╔═╡ 33e6b93d-ad9d-4e58-93f9-3e26b4b579d9
if typeof(A) == Array{Float64, 2}
	vectorfield(A, if n normalize else identity end, scale)
end

# ╔═╡ 3e7ba883-2735-4edd-8f30-1cbf59c265f0
vectorfield(J(ϵ), if nJ normalize else identity end, scaleJ)

# ╔═╡ Cell order:
# ╟─c44dff8f-8e7a-4904-9a74-2cbf39250a03
# ╟─8dc480f7-d8be-4826-bc0b-187ec47cf14d
# ╟─2a237878-db73-4f29-87ea-cbcc49b642d0
# ╟─e75765e3-85fa-4eb1-abb0-7d740ec091cb
# ╟─1c12a305-572d-44e5-8f25-7b0ecfc1313b
# ╟─bfac1751-6987-431b-967e-3305df530308
# ╟─86f58406-08c7-4adb-95ae-524c0d4666e7
# ╟─33e6b93d-ad9d-4e58-93f9-3e26b4b579d9
# ╟─08a1e247-4ba9-4693-bf11-06cc519320fc
# ╟─69451070-2f2a-4052-8b97-40625819db00
# ╟─5479c733-f7c7-4ec7-8b4c-baff5fba9561
# ╟─67adecf7-21e2-43aa-9dea-f881259ed939
# ╟─2b4e0633-025a-4ec9-8a01-24e11e1e8fd5
# ╟─1d59a29e-a930-481f-8093-3a232a1b6962
# ╟─54d73ea3-86fe-4f7b-958a-6090e571bea2
# ╟─3e7ba883-2735-4edd-8f30-1cbf59c265f0
# ╟─f5a77a25-cc64-4d08-93fd-42beb3007cb4
# ╟─fe63c64e-fef5-4abc-9636-d1413d40fc1b
# ╟─a3e0aee0-b09c-4980-85b8-066cf9881757
# ╟─d21e963c-a86a-11eb-3082-41cb8d8b1c17
# ╟─8f29d603-c929-4004-a1f9-afd59e66b72f
# ╟─c9edac8f-8b56-4402-9315-9421f2bfbc7f
# ╟─8edbf3a3-6f75-4eff-887f-07005b860f43
