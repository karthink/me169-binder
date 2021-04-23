### A Pluto.jl notebook ###
# v0.14.2

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

# ╔═╡ 439697ba-a406-11eb-1588-ebb4008dad62
begin
	using Plots, PlutoUI, LinearAlgebra, Printf
end

# ╔═╡ e6296659-b446-4c0b-9a91-123c98ed3862
begin
	Base.@kwdef mutable struct Firefly
		dt = 0.02
		α = pi/3.0
		θ = 0
		A = 0.9
		
		Ω = 2.8
		ω = 1.0
	end
	
	function step!(f::Firefly)
		f.α += f.Ω * f.dt
		f.θ += (f.ω + f.A * sin(f.α - f.θ)) * f.dt
		
# 		f.α = f.α % (2*pi)
# 		f.θ = f.θ % (2*pi)
	end
	
	angles(f::Firefly) = [ (cos(f.θ), sin(f.θ)), (cos(f.α), sin(f.α)) ]
end

# ╔═╡ 0b1976a5-ed63-4ad7-bb97-2ade4f8ab3d0
function phase_animate(f::Firefly)
	@gif for i in 1:450
		plt = plot(cos, sin, 0, 2*pi,
			    axis=false, ticks=false, grid=false,
				aspect_ratio=1.0, label=false)
		scatter!(plt, angles(f), 
			     color=[:red, :blue], 
			     markersize=7.0, label=false)
		ϕ_str = @sprintf("%4f", f.α - f.θ)
		annotate!(plt, [(1,1,text("Firefly",:top,:right,:red)),
				        (1,1,text("Flashlight",:bottom,:right,:blue)),
				        (0.2*cos(0.5*(f.α + f.θ)),
					     0.2*sin(0.5*(f.α + f.θ)),
					     text("ϕ"; fontsize=15)),
				        (-1,1,text("ϕ = $(ϕ_str)";fontsize=12))])
		plot!(plt, t -> t*cos(f.α), t-> t*sin(f.α), 0, 1,
			  linestyle=:solid, linecolor=:blue, label=false)
		plot!(plt, t -> t*cos(f.θ), t-> t*sin(f.θ), 0, 1,
			  linestyle=:dash, linecolor=:red, label=false)
		if abs((f.α % (2*pi)) - pi/2) < 0.05
			plot!(plt, cos, sin, 0, 2*pi,
				  fill = true, color = :blue,
				  alpha = 0.35, label=false)
		end
		if abs((f.θ % (2*pi)) - pi/2) < 0.05
			plot!(plt, cos, sin, 0, 2*pi,
				fill = true, color = :red, 
				alpha = 0.35, label=false)
		end
		step!(f)
	end
end

# ╔═╡ 1f78957f-7c62-495d-8ae8-91a8364775fe
function fplot(mu)
	fplot = plot(;label=false, framestyle=:origin)
	plot!(fplot, t -> mu - 0.9*sin(t), xlims = (-2*pi, 4*pi), color=:blue, 	                     linewidth=2.5, label=false, 
		  ylabel = "$(@sprintf("%0.2f",mu)) - sin(ϕ)", xlabel="ϕ",
		  xticks = range(-6,12, step=2)
		  	)
	
end

# ╔═╡ dd7ba46a-ba61-4193-8ccf-b73d39ce520f
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))

# ╔═╡ bb317040-c683-4674-9f08-855c327897da
md"""
## Flows on a Circle

By flowing in one direction, the state can eventually return to its starting place - periodic motion in a 1D system!

### Example: Uniform oscillator
```math
\dot{\theta} = \Omega
```
"""

# ╔═╡ 4b5fd975-61e7-481f-aab7-f4c78c724a7d
md"""
### Is the state space a circle?

What is wrong with the system
```math
\dot{\theta} = \theta,\quad \theta \in \mathbb{R}
```
on a circle?

**Dynamics on a circle: $\dot{\theta} = f(\theta)$ needs to have a special property:**
"""

# ╔═╡ c39dfeea-8764-41c7-9cdc-a1c874f85b64
hint(md"""
	$\dot{\theta} = f(\theta),\quad f(\theta + 2 \pi) = f(\theta)$
	Otherwise $\dot{\theta}$ is not well-defined!
	""")

# ╔═╡ 4ff1a8fb-9a84-4bfe-bf95-10fb3587ff52
let
	straight = plot( identity, xlims = (0, 4*pi), 
	xlabel = "θ", ylabel = "f(θ)", lw=2.5,
	label = "f(θ)")
	vline!(straight, [1.5, 1.5 + 2*pi], ls=:dash, lc=:red, label=false, lw=2.5)
	plot!(straight, [(t,10) for t in range(1.5,1.5+2*pi,step=0.01)], arrows=true, label=false, lw=1.5, lc=:red)
	annotate!(straight, [(1.5+pi,10,text("2π",:top,:red))])
end

# ╔═╡ e5182409-101e-4b6a-8282-3faf443fdc7e
md"""
## Fireflies
"""

# ╔═╡ fd5fd2c9-c654-453e-bf01-6280035e58c0
html"""
<iframe width="560" height="315" src="https://www.youtube.com/embed/EnwVVE-EGVw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ d71a14cf-6c1d-4674-81ed-254e85554365
md"""
Our model is (_see Strogatz, section 4.5_)
```math
\begin{align*}
\dot{\alpha} &= \Omega\\
\dot{\theta} &= \omega + A \sin(\alpha - \omega)\\
\phi &:= \alpha - \theta\\
\end{align*}
```
Which we rewrite in terms of $\phi$ as
```math
\begin{align*}
\dot{\alpha - \theta} &= \Omega - \omega - A \sin(\alpha - \omega),\ \text{or}\\
\implies \dot{\phi} &= (\Omega - \omega) - A \sin(\phi)
\end{align*}
```

When nondimensionalized, this gives us
```math
\frac{d \phi}{d \tau} = \mu - sin(\phi)
```
where $\tau = A t$, $\mu = \frac{\Omega - \omega}{A}$
"""

# ╔═╡ 960fefd7-79c8-4b62-a9f9-0bef0954268e
md"""
### Case 1: Identical flashing frequencies, $\Omega = \omega$"""

# ╔═╡ 2c711db1-6e30-4d72-981d-93b5c9977f42
fplot(0)

# ╔═╡ f2da5d6a-19ad-48ca-ae98-fe84207cf544
md"**Click this checkbox to see what the dynamics look like!**"

# ╔═╡ 08d6901e-4468-4f81-b6ae-22ee5635d779
md"""Case: $\Omega = \omega = 1$ $(@bind c11 html"<input type=checkbox >")"""

# ╔═╡ e786a03a-ecfa-4add-9e4c-fad6755c43aa
if c11
	phase_animate(Firefly(Ω=1.0, ω=1.0))
end

# ╔═╡ 3c049ac7-4255-4ec6-872d-697e44ceab4a
md"""
### Case 2: $\Omega > \omega,\ \lvert\Omega - \omega\rvert \le A$
"""

# ╔═╡ 32831825-f8de-4801-b7ec-7c7b2fabdf2b
fplot((1.5-1))

# ╔═╡ 04667541-a8d2-4d86-828d-2dbf8fde5301
md"""Case: $\Omega = 1.5,\ \omega = 1$ $(@bind c151 html"<input type=checkbox >")"""

# ╔═╡ aca698ef-a475-454c-a9f2-1a6bbe984522
if c151
	phase_animate(Firefly(Ω=1.5, ω=1))
end

# ╔═╡ 583dbde8-7808-46bc-93b3-cae9562e1469
fplot((1.9-1))

# ╔═╡ 57192cb9-042a-47d1-9273-584e3b123952
md"""Case: $\Omega = 1.9,\ \omega = 1$ $(@bind c191 html"<input type=checkbox >")"""

# ╔═╡ ba010d9d-2b9a-4a04-8dff-287f082088e7
if c191
	phase_animate(Firefly(Ω=1.9, ω=1))
end

# ╔═╡ 7a6e3ba5-630d-4c57-a836-8b1cbabad181
md"""
### Case 2: $\Omega > \omega,\ \lvert\Omega - \omega\rvert > A$
"""

# ╔═╡ 8be21c48-34a3-4d52-8d58-dfc6c9c64aff
fplot((2.0-1))

# ╔═╡ fc16376f-9f44-47a8-bddd-0b695eb1a239
md"""Case: $\Omega = 2.0,\ \omega = 1$ $(@bind c251 html"<input type=checkbox >")"""

# ╔═╡ 8aab9b1a-38a6-40d6-8ab1-8e59daf3634b
if c251
	phase_animate(Firefly(Ω=2.0, ω=1, θ=0, α=pi/2))
end

# ╔═╡ Cell order:
# ╟─439697ba-a406-11eb-1588-ebb4008dad62
# ╟─e6296659-b446-4c0b-9a91-123c98ed3862
# ╟─0b1976a5-ed63-4ad7-bb97-2ade4f8ab3d0
# ╟─1f78957f-7c62-495d-8ae8-91a8364775fe
# ╟─dd7ba46a-ba61-4193-8ccf-b73d39ce520f
# ╟─bb317040-c683-4674-9f08-855c327897da
# ╟─4b5fd975-61e7-481f-aab7-f4c78c724a7d
# ╟─c39dfeea-8764-41c7-9cdc-a1c874f85b64
# ╟─4ff1a8fb-9a84-4bfe-bf95-10fb3587ff52
# ╟─e5182409-101e-4b6a-8282-3faf443fdc7e
# ╟─fd5fd2c9-c654-453e-bf01-6280035e58c0
# ╟─d71a14cf-6c1d-4674-81ed-254e85554365
# ╟─960fefd7-79c8-4b62-a9f9-0bef0954268e
# ╠═2c711db1-6e30-4d72-981d-93b5c9977f42
# ╟─f2da5d6a-19ad-48ca-ae98-fe84207cf544
# ╟─08d6901e-4468-4f81-b6ae-22ee5635d779
# ╟─e786a03a-ecfa-4add-9e4c-fad6755c43aa
# ╟─3c049ac7-4255-4ec6-872d-697e44ceab4a
# ╠═32831825-f8de-4801-b7ec-7c7b2fabdf2b
# ╟─04667541-a8d2-4d86-828d-2dbf8fde5301
# ╟─aca698ef-a475-454c-a9f2-1a6bbe984522
# ╠═583dbde8-7808-46bc-93b3-cae9562e1469
# ╟─57192cb9-042a-47d1-9273-584e3b123952
# ╟─ba010d9d-2b9a-4a04-8dff-287f082088e7
# ╟─7a6e3ba5-630d-4c57-a836-8b1cbabad181
# ╠═8be21c48-34a3-4d52-8d58-dfc6c9c64aff
# ╟─fc16376f-9f44-47a8-bddd-0b695eb1a239
# ╟─8aab9b1a-38a6-40d6-8ab1-8e59daf3634b
