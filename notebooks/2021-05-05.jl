### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ cfdc39ea-ade9-11eb-2036-2fc71efeb18d
begin
	using Plots
	using Pluto
	using PlutoUI
	using LinearAlgebra, DifferentialEquations
end

# ╔═╡ 86790045-1149-41f6-ad26-88bbd57bbab3
begin
function vectorfield(f; norm=identity, scale=0.2, params=nothing, xlims=(-2.0, 2.0), ylims=(-2.0, 2.0))
	plt = plot(framestyle=:origin, aspect_ratio=1.0, xlims = xlims, ylims = ylims)
	# plotvec(plt,eigvecs(A))
	gridnum = 25

	x = range(xlims[1], stop=xlims[2], length=gridnum)
	y = range(ylims[1], stop=ylims[2], length=gridnum)
	f_marker = (t -> scale .* t) ∘ norm ∘  ((x,y) -> f((x,y), params, nothing))
	quiver!(repeat(x,gridnum), vec(repeat(y',gridnum)), quiver=f_marker, c=:blue, arrowscale=0.1, alpha=0.5, lc="#1d1f21")
	
	return plt
end

end

# ╔═╡ 948c2f8e-a0f2-4521-812c-a02758e0ecee
begin
	function trajectories(f;tfinal=3.0, 
		           params=nothing,
				   xstart=range(-2.0,2.0,length=5),
			       ystart=range(-2.0,2.0,length=5),
		           xlims=(-2.0,2.0),
		           ylims=(-2.0,2.0))
		plt = plot(framestyle=:origin, aspect_ratio=1.0, xlims = xlims, ylims = ylims, xticks=false, yticks=false, label=false)
    
      for x0 in xstart
		for y0 in ystart
			prob = ODEProblem(f, [x0 ; y0], (0, tfinal), params)
        	sol = solve(prob)
        	plot!(plt, sol, vars=(1,2), lw = 1.5, label=false, arrow=true, color=:black, aspect_ratio=1.0, xlims = xlims, ylims = ylims)
      end
		end
      return plt
	end
end

# ╔═╡ 7a230eac-a55a-42d9-a60f-ff7f32c5892a
md"""
### Vector field & Trajectory calculator

To plot the vector field $f(x,y)$, modify this definition of f:
"""

# ╔═╡ 52d2a5ee-c1d7-4027-86d8-b55be45ad6ce
f((x,y), p, t) = [ -y; +x ]

# ╔═╡ 7be48595-b380-4fe8-a040-6273bfdf9c49
md"You can use params as a parameter of the equation. For instance, with parameters $r$ and $s$ in
```math
f(x,y) = \begin{bmatrix} - r y \\ s x \end{bmatrix}
```
you can write `f((x,y), (r,s), t) = [ -r*y; s*x ]`"

# ╔═╡ a9c6d309-83fc-41b6-a8eb-8beeda11b598
md"Set your system parameters (if you have any) here:"

# ╔═╡ a87c4661-2986-4d94-b9a0-458f00b0c1e2
r,s = 1,1

# ╔═╡ 2d5f2740-80b5-41ca-afaa-886a3a202ac2
md"""
#### Vector field calculator

Modify the following function call to suit your needs:
"""

# ╔═╡ 939ff4b1-4678-4156-a728-cd11708f14d0
vectorfield(f, 
	        params=(r,s),      # If your system has no params use params=nothing
			xlims=(-1.0, 1.0), # x-limits of plot
	        ylims=(-1.0, 1.0), # y-limits of plot
			norm=identity,     # Use norm=normalize to make all arrows the same length
	       	scale=0.2)         # Use to make arrows smaller/larger

# ╔═╡ cb9b459e-b67a-4ffc-92a3-8119350a531c
md"""
#### Trajectory calculator

Modify the following function call to suit your needs:
"""

# ╔═╡ 83cc9479-f099-49ab-9dbe-99c0399465f3
trajectories(f,
	tfinal=3.0,                      # Simulate until time tfinal
	params=(r,s),                    # If your system has no params use params=nothing
	xlims=(-1.0, 1.0),               # x-limits of plot
	ylims=(-1.0, 1.0),               # y-limits of plot
	xstart=range(-1.0,1.0,step=0.4), # trajectories starting from x,y = xstart, ystart
	ystart=range(-1.0,1.0,step=0.4)) # trajectories starting from x,y = ystart, ystart

# ╔═╡ 733685a6-749a-4f72-a07a-c5ae2abbac20
md"""
### Plots for discussion session, week 6
"""

# ╔═╡ cbd78270-795c-4c36-9996-0355c30868e2
md"""
##### Strogatz 6.3.10

```math
f(x) = \begin{bmatrix} x y \\ x^2 - y \end{bmatrix}
```
"""

# ╔═╡ cd7b0468-60fc-4075-b9f9-a91f09e9d57e
md"##### Vector field and nullclines:"

# ╔═╡ d2cc4fea-cc01-417b-975e-813d1ae04f28
let
	f((x,y), p, t) = [x*y; x^2 - y]
	plt = vectorfield(f; norm=identity, scale=0.2)
	plot!(plt, t -> t^2, label=false, lc=:blue, lw=2.0)
	vline!(plt, [0], label=false, lc=:red, lw=2.0)
	hline!(plt, [0], label=false, lc=:red, lw=2.0)
end


# ╔═╡ 3a1d61fd-5295-4906-9595-0d3848f6973a
md"##### Trajectories:"

# ╔═╡ c46bdfbd-f44f-4a4d-b1fd-fcd830b57d0f
let
	f((x,y), p, t) = [x*y; x^2 - y]
	plt = trajectories(f, tfinal=3.0, params=nothing, xstart=range(-1,1,length=7), ystart=range(-1,1,length=7))
	plot!(plt, xlims=(-1,1), ylims=(-1,1), xticks=true, yticks=true)
end

# ╔═╡ c8a8f82f-a408-4fa8-bd57-1009af7f2ba0
md"""
##### Strogatz 6.4.11

```math
\begin{align}
\dot{l} &= \alpha l c \\
\dot{r} &= \alpha r c \\
\dot{c} &= - \alpha l c - \alpha r c
\end{align}
```
"""

# ╔═╡ 53a8e6b5-8eb7-4801-9069-f9f3292aee56
pol((x,y),α,t) = [ α*x*(1-x-y); α*y*(1-x-y) ]

# ╔═╡ d1098f1b-00aa-49a8-9abc-c8ec8ae53b20
md"##### Vector Field and nullclines"

# ╔═╡ d8489248-0c2b-4ed5-a408-d4cfb3f5a4f2
let
	plt = vectorfield(pol, params=-1.0, norm=identity, scale=0.2 )
	plot!(plt, xlims=(0,1.5), ylims=(0,1.5))
	plot!(plt, t -> 1 - t, xlims=(0,1.5), ylims=(0,1.5), lc=:purple, lw=2.0)
	hline!(plt, [0], lc=:red, lw=2.0)
	vline!(plt, [0], lc=:blue, lw=2.0)
end

# ╔═╡ 6b03901b-3762-40a1-a1c5-04531781a39b
md"##### Trajectories"

# ╔═╡ 14259ca8-cde9-41e9-a3e7-7940595b079b
let
	plt2 = trajectories(pol, params=1.0, tfinal=3.0,
	                    xstart=range(0,0.7,length=5),
		                ystart=range(0,0.7,length=5))
	plot!(plt2, xlims=(0,1.5), ylims=(0,1.5))
end

# ╔═╡ 4d78ceb2-c9e7-4ecd-adf6-3ac3f17dfc73
md"""
##### Strogatz 6.7.2

```math
\begin{align}
\ddot{\theta} = γ - \sin{\theta}
\end{align}
```
"""

# ╔═╡ 5e17f5fc-0d2c-4a8c-9258-62bb8cb4b6ce
γ = 0.5

# ╔═╡ 0a010364-b9e2-494d-a3c7-48ae4afb28fd
md"##### Vector Field and nullclines"

# ╔═╡ d11dc629-674d-4a01-aa4b-4d18202b1445
let
	pend((x,y), γ, t) = [ y, γ - sin(x) ]
	lims = (-5.0,5.0)
	plt = vectorfield(pend, params=γ, xlims=lims, ylims=lims)
	vline!(plt, [γ], lc=:teal, lw=2.0, xticks=true, xlims=(-3.0,5.0))
	hline!(plt, [0], lc=:red,  lw=2.0)
end

# ╔═╡ 22b6162d-ac66-4e35-bc96-4e3b3e452562
md"##### Trajectories"

# ╔═╡ 574d50d5-d432-463a-9e96-81e8df666e83
let
	pend((x,y), γ, t) = [ y, γ - sin(x) ]
	xr = (-1.0, 20.0)
	yr = (-4.0, 10.0)
	plt2 = trajectories(pend, tfinal=6.0, params=γ, xlims=xr, ylims=yr, xstart=range(-1.0,10.0,length=12), ystart=range(-1.0,5.0,length=5))
	plot!(plt2, xlims=xr, ylims=yr)
end

# ╔═╡ Cell order:
# ╟─cfdc39ea-ade9-11eb-2036-2fc71efeb18d
# ╟─86790045-1149-41f6-ad26-88bbd57bbab3
# ╟─948c2f8e-a0f2-4521-812c-a02758e0ecee
# ╟─7a230eac-a55a-42d9-a60f-ff7f32c5892a
# ╠═52d2a5ee-c1d7-4027-86d8-b55be45ad6ce
# ╟─7be48595-b380-4fe8-a040-6273bfdf9c49
# ╟─a9c6d309-83fc-41b6-a8eb-8beeda11b598
# ╠═a87c4661-2986-4d94-b9a0-458f00b0c1e2
# ╟─2d5f2740-80b5-41ca-afaa-886a3a202ac2
# ╠═939ff4b1-4678-4156-a728-cd11708f14d0
# ╟─cb9b459e-b67a-4ffc-92a3-8119350a531c
# ╠═83cc9479-f099-49ab-9dbe-99c0399465f3
# ╟─733685a6-749a-4f72-a07a-c5ae2abbac20
# ╟─cbd78270-795c-4c36-9996-0355c30868e2
# ╟─cd7b0468-60fc-4075-b9f9-a91f09e9d57e
# ╟─d2cc4fea-cc01-417b-975e-813d1ae04f28
# ╟─3a1d61fd-5295-4906-9595-0d3848f6973a
# ╟─c46bdfbd-f44f-4a4d-b1fd-fcd830b57d0f
# ╟─c8a8f82f-a408-4fa8-bd57-1009af7f2ba0
# ╠═53a8e6b5-8eb7-4801-9069-f9f3292aee56
# ╟─d1098f1b-00aa-49a8-9abc-c8ec8ae53b20
# ╟─d8489248-0c2b-4ed5-a408-d4cfb3f5a4f2
# ╟─6b03901b-3762-40a1-a1c5-04531781a39b
# ╟─14259ca8-cde9-41e9-a3e7-7940595b079b
# ╟─4d78ceb2-c9e7-4ecd-adf6-3ac3f17dfc73
# ╠═5e17f5fc-0d2c-4a8c-9258-62bb8cb4b6ce
# ╟─0a010364-b9e2-494d-a3c7-48ae4afb28fd
# ╟─d11dc629-674d-4a01-aa4b-4d18202b1445
# ╟─22b6162d-ac66-4e35-bc96-4e3b3e452562
# ╟─574d50d5-d432-463a-9e96-81e8df666e83
