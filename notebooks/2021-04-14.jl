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

# ╔═╡ 9d60a7b0-9d92-11eb-1281-1bd6c9aacf25
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	using Plots, PlutoUI
end

# ╔═╡ b9b651a4-f540-4969-a5c4-3da849492db0
md"""
Recall from the discussion session that we reduced the overdamped bead on a hoop problem to the ODE
```math
\epsilon \ddot{\phi} + \dot{\phi} = f(\phi) = \sin(\phi) \left( \gamma \cos(\phi) - 1 \right)
```

We are looking at this system in the limit $\epsilon \to 0$, so that we may approximate it by the first order system
```math
\dot{\phi} = f(\phi) = \sin(\phi) \left( \gamma \cos(\phi) - 1 \right)
```

But approximating a second order system with a first order one is tricky business: Our reduced system is now overdetermined by specifying initial conditions to the original equation. Given arbitrary $(\phi(0), \dot{\phi}(0))$, _we cannot find a solution to the first order equation_. 

Understanding the relationship between the solutions of these two equations is the purview of singular perturbation theory, which we examine in the following section.
"""

# ╔═╡ 76711d3e-efdf-403b-b8d5-e114e42c0adf
md"""
## Singular Perturbation Theory

To understand singular perturbation theory applied to the overdamped bead on the rotating hoop problem (Section 3.5 in Strogatz), we're going to start by looking at a simpler system we already fully understand: The RC circuit from week 1.

### Capacitor charging revisited

You may recall that the ODE that describes the charge dynamics is
```math
\dot{Q} = \frac{V}{R} - \frac{Q}{R C},\quad Q(0) = Q_0
```

We want to find conditions under which we can neglect the highest order derivative. Since this is a first order system, a "lower order" approximation would be just an algebraic equation!

Specifically, it would be the algebraic system
```math
0 = \frac{V}{R} - \frac{Q(t)}{R C}
```
which has the constant solution $Q(t) = C V$. You may recognize this as the steady state solution of the differential equation starting from any initial condition.

So the salient questions are:

1. Under what conditions can we approximate the solution of the $1^{st}$ order ODE by the solution of an algebraic equation?
2. How do we reconcile starting with an arbitrary initial $Q_0$ (which the original ODE allows for) with the solution of the algebraic equation (which does not)?

These questions are the purview of singular perturbation theory. To answer them, we will first rescale the ODE to be nondimensional. This makes it easier to compare the magnitudes of the three terms in the equation.

Define $q := \frac{Q}{CV}$. Then the equation becomes
```math
R C \dot{q} = 1 - q,\quad q(0) = \frac{Q₀}{CV}
```

Where all three terms are now nondimensional. As we've seen before, it produces solutions like the following:
"""

# ╔═╡ 67fde8e0-3e5e-4a24-acf0-91bb8e7b82c8
md"""
*Move me*: q₀ = $(@bind q₀ Slider(0.0:0.05:2.0; default=0.25, show_value=true))
"""

# ╔═╡ d08e065a-acd7-4f73-bd3b-5fb1958233e6
let
	Q(t) = 1 - (1 - q₀) * exp(- t)
	plot(Q, xlims = (0.0, 10.0), ylims = (0, 2.5), 
		 linewidth=2.0, ticklabels=false,
	     xlabel="t", ylabel="q", label="q(t)")
	hline!([1.0], linestyle=:dash, linecolor=:green, alpha=0.7, label=false)
	annotate!([(0.0, q₀, text("q₀ = $(q₀)", 13, :left, :top, :blue))])
end

# ╔═╡ 758a737c-f0fe-4e74-a2ed-c8067999b93d
md"""
We would like to compare the relative magnitudes of the three terms in the equation to determine which one(s) we can safely neglect. Of the three, the RHS terms are $O(1)$ by construction. But the derivative term ($RC \dot{q}$) can have varying magnitudes at different times.

The piece of the puzzle that's missing is the *time scale* at which we examine the solution.

Let's rescale time as $\tau = \frac{t}{T}$, where $T$ is a parameter we can tune. From the chain rule, 
```math
  \frac{d}{d t} = \frac{1}{T} \frac{d}{d \tau}
```
and our nondimensional equation becomes
```math
\frac{R C}{T} \frac{d q}{d \tau} = 1 - q
```

----
**What does it mean to rescale time in a differential equation?** 

Each unit of nondimensional time ($\tau$) corresponds to $T$ units of actual time $t$, so you can think of rescaling $t$ as zooming into (or out of) the graph of a function $f(t)$. As you change the zoom level, the function (and the ODE it obeys) can look remarkably different over the same range of $\tau$. Observe the change in the rescaled function $f(\tau)$ corresponding to the original $f(t) = \sin(t)/t$ as we vary the "zoom level" $S$ in $\tau := t / S$: """

# ╔═╡ fccf0efa-a831-4627-819c-711b0be6e762
@bind spow Slider(-6:6, default = -4)

# ╔═╡ 7a25ea0a-ee74-4959-9518-eddf9902534b
begin
	S = 2.0^spow
	md" $S$ = $(S)"
end

# ╔═╡ 008e3bd5-522d-407f-a2d0-12be82069056
begin
	function f(t::Real)
		if t ≠ 0
			sin(t * S)/(t * S)
		else
			1.0
		end
	end
	plot(τ -> f(τ), xlims = (0,40), ylims = (-1.3, 1.3),
			label = false, linewidth = 2.0, alpha = 0.9,
	     xlabel = "τ", ylabel = "f(τ)")
end

# ╔═╡ 20b1f5cf-4e4c-4a0d-abe5-6395f9503dd9
md"""
Alternatively, you can think of it as squeezing or stretching the graph of $f(\tau)$ in the $\tau$-direction.

----
"""

# ╔═╡ 9b87bb42-9f10-462f-abcf-f1de26963086
md" Now study the effect of varying the parameter $\frac{RC}{T}$ on the solution $q(\tau)$ of
```math
\frac{R C}{T} \frac{d q}{d \tau} = 1 - q
```
"

# ╔═╡ f451c7cd-6a87-447e-99b1-99d4542eba46
begin
	@bind pow Slider(-6:10; show_value=false, default=2)
end

# ╔═╡ 1f456fdf-9519-4737-8926-c904b648b365
begin
	RCbyT = 2.0^pow
	md"``\frac{RC}{T} = `` $(RCbyT)"
end

# ╔═╡ 2bad47ce-c31e-471b-91bd-7e248e907e56
let
	q(t) = 1 - (1 - q₀)*exp(-t / (RCbyT))
	plot(q, xlims = (0, 40.0), ylims = (0, 2.5),
		 linewidth = 2.0, label = "q(τ)",
		 xlabel = "τ", ylabel="q(τ)")
end

# ╔═╡ 6474de40-b857-4aa9-8a7c-e6f309e79123
md"""
1. At which extreme of $\frac{RC}{T}$ are we justified in neglecting the derivative term (thus reducing the system to an algebraic equation)? What does this imply for the time constant of the circuit?
2. In this limit, how does the system reconcile starting from an initial condition that is different from $q = 1$?

#### Case $\frac{RC}{T} \ll 1$
We can gain more insight into the fast transient at $\tau = 0$ when $\frac{RC}{T} \to 0$, or for large $T$. The nondimensional equation at $\tau = 0$ is
```math
\frac{d q}{d \tau}_{\tau = 0} = \frac{T}{R C} (1 - q_0)
```
If $q_0 < 1$, $\frac{d q}{d \tau}$ is large and positive, so the trajectory shoots up until the right side is of $O(1)$, or until $(1 - q) \sim \frac{R C}{T}$. The behavior for $q_0 > 1$ is similar. Try changing $\frac{RC}{T}$ above:
"""

# ╔═╡ ce6d606d-77f1-4ed4-b957-204ade5acecb
p = let
	q(q0, t) = 1 - (1 - q0)*exp(-t / (RCbyT))
	p = plot(xlims = (0,40.0), ylims = (0.0,2.0),
	         xlabel= "τ", ylabel="q(τ)")
	for k in range(0,1.75,step=0.21)		
	plot!(p, t -> q(k,t), label=false, alpha=0.70, linewidth=2.0)
	end
	p
end

# ╔═╡ 927f272e-3eb6-49e1-a275-4cdf604eace5
md"""
Thus for small $\frac{RC}{T}$, all solutions slam into the line $1 - q = 0$, and this happens in time $T \sim \frac{RC}{1-q_0} \sim RC$.

#### Case $\frac{RC}{T} \gg 1$

In this case the derivative at $\tau = 0$
```math
\frac{d q}{d \tau}_{\tau = 0} = \frac{T}{RC} (1 - q_0)
```
is essentially zero, since $\frac{T}{RC} \to 0$ and $(1 - q_0) \sim O(1)$, and the solutions barely change. (Try moving the $\frac{RC}{T}$ slider above all the way to the right.)
"""

# ╔═╡ 8ae2a65b-db49-43b2-a3ca-0557fb3cea36
p

# ╔═╡ 5e46a9b6-03c8-4d43-a29b-cb993b688faa
md"""
#### Case $\frac{RC}{T} \sim 1$
This is the interesting case, now all terms in the nondimensional equation are $O(1)$:
```math
\frac{RC}{T} \frac{d q}{d \tau} = \frac{d q}{d \tau} = 1 - q
```
At this timescale we cannot neglect any of the terms, and we would be ill-advised to approximate the solution of the ODE by dropping the derivative entirely. $T = RC$ is the timescale over which the transient behavior occurs, and as you might have expected, this is exactly the time constant of an RC circuit.
"""

# ╔═╡ 319a0160-b574-4528-99ec-15ebe034d419
md"""
#### Singular perturbation of the capacitor charging problem: Summary

Here are the questions we posed at the beginning of the section:

**Q: Under what conditions can we approximate the solution of the  order ODE by the solution of an algebraic equation?**

We can do this on when our time scale of interest $T$ such that $T \gg RC$, _i.e._ our time scale of interest is much larger than the time constant of the circuit.

**Q: How do we reconcile starting with an arbitrary initial  (which the original ODE allows for) with the solution of the algebraic equation (which does not)?**

On this long time scale, there is a fast initial transient (corresponding to time $\sim RC$) during which all initial conditions essentially collapse to the solution of the algebraic equation. Of course, this solution is the fixed point of the ODE.
"""

# ╔═╡ fa0beecd-3519-472f-afa4-6c2590c8df09
md"""
### Singular perturbation of the overdamped bead on a rotating hoop

The RC circuit problem was, in a sense, trivial to apply singular perturbations to, since the reduced equation simply gives us the fixed point of the system. But we can apply the same idea to the overdamped bead problem:
```math
m r \ddot{\phi} = - b \dot{\phi} + m g f(\phi)
```
where $f(\phi) = \sin{\phi} \left( \frac{r \omega^2}{g} \cos{\phi} - 1 \right) =: \sin{\phi} ( \gamma \cos{\phi} - 1 )$

In analyzing this ODE as a first order system, we neglected the inertia term with the assumption that we are justified in doing so. Here we will find the conditions (including the timescale) under which this assumption is fair.
"""

# ╔═╡ 174a248d-1033-44f7-83d8-d491384bdd4d


# ╔═╡ 18fe0964-44ff-4456-a275-a9dd68a381d8
md"### Helper code"

# ╔═╡ Cell order:
# ╟─b9b651a4-f540-4969-a5c4-3da849492db0
# ╟─76711d3e-efdf-403b-b8d5-e114e42c0adf
# ╟─67fde8e0-3e5e-4a24-acf0-91bb8e7b82c8
# ╟─d08e065a-acd7-4f73-bd3b-5fb1958233e6
# ╟─758a737c-f0fe-4e74-a2ed-c8067999b93d
# ╟─fccf0efa-a831-4627-819c-711b0be6e762
# ╟─7a25ea0a-ee74-4959-9518-eddf9902534b
# ╟─008e3bd5-522d-407f-a2d0-12be82069056
# ╟─20b1f5cf-4e4c-4a0d-abe5-6395f9503dd9
# ╟─9b87bb42-9f10-462f-abcf-f1de26963086
# ╟─f451c7cd-6a87-447e-99b1-99d4542eba46
# ╟─1f456fdf-9519-4737-8926-c904b648b365
# ╟─2bad47ce-c31e-471b-91bd-7e248e907e56
# ╟─6474de40-b857-4aa9-8a7c-e6f309e79123
# ╟─ce6d606d-77f1-4ed4-b957-204ade5acecb
# ╟─927f272e-3eb6-49e1-a275-4cdf604eace5
# ╟─8ae2a65b-db49-43b2-a3ca-0557fb3cea36
# ╟─5e46a9b6-03c8-4d43-a29b-cb993b688faa
# ╟─319a0160-b574-4528-99ec-15ebe034d419
# ╟─fa0beecd-3519-472f-afa4-6c2590c8df09
# ╠═174a248d-1033-44f7-83d8-d491384bdd4d
# ╟─18fe0964-44ff-4456-a275-a9dd68a381d8
# ╟─9d60a7b0-9d92-11eb-1281-1bd6c9aacf25
