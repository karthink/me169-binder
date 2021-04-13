### A Pluto.jl notebook ###
# v0.12.21

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

# ╔═╡ 7b18828e-984e-11eb-0f5a-ff870f4b3f7e
begin
	using PlutoUI, Plots
end

# ╔═╡ e43df1dc-9735-11eb-2361-8dbf881cd29d
md"## The Forward Euler Method"

# ╔═╡ afac1ac4-97ef-11eb-1521-e36237bb6988
md"""
Applied to
```math
{\huge \dot{x} = f(x)}
```
with ``x(0) = x_0``, and ``t_0 \le t \le t_1``.
"""

# ╔═╡ 6dd49d18-982d-11eb-3a78-8bb79e6cd97a
md"""
The Forward Euler method is the most straightforward explicit method of integrating this system numerically. We find approximate values of ``x(t)`` at evenly spaced intervals ``n \Delta t`` by using successive linear approximations of the function. 

Denoting the numerical solution by ``x[\cdot]``, the first step approximates ``x(\Delta t)`` by ``x[\Delta t]`` as:
```math
x[1 \Delta t] = x_0 + f(x_0) \Delta t \approx x_0 + f(x_0) \Delta t + O(\Delta t^2)  = x(\Delta t)
```
At each following step, we discard the (unknown) ``O(\Delta t^2)`` terms:

```math
\begin{align*}
x[1 \Delta t] &= x_0 + f(x_0) \Delta t\\
x[2 \Delta t] &= x[1 \Delta t] + f(x[1 \Delta t]) \Delta t\\
\vdots &= \vdots\\
x[(n+1) \Delta t] &= x[n \Delta t] + f(x[n \Delta t]) \Delta t\\
\end{align*}
```

In the below examples, try playing with the initial condition and step size ``\Delta t`` to get an idea of the accuracy of this method.
"""

# ╔═╡ 75050150-97f0-11eb-10dd-45fff81be40b
md"""
---
### Capacitor charging: $\dot{x} = f(x) := 1 - x$
"""

# ╔═╡ 18a72e4e-9738-11eb-30c2-a5700e9ac19f
@bind pow_capacitor Slider(-10.0:0.2:2.0, default=0.6, show_value=true)

# ╔═╡ c190a43c-9739-11eb-071d-a721c779ba74
begin
	Δt_capacitor = 2.0^pow_capacitor
	md"""
	#### Step size:
	
	The step size ``{\large \Delta t_{capacitor} =}`` $Δt_capacitor.
	
	Change it using the slider below!
	"""
end

# ╔═╡ abdbbe84-9814-11eb-36a8-2b769e0eae65
md"""
#### Starting from
``{\large x_0}`` = $(@bind x₀_capacitor Slider(0:0.01:2, default=0.25,show_value=true))
"""

# ╔═╡ 3c63472a-982f-11eb-361a-355745c9d948
md"The below plot shows the error, which is the difference between the true solution - and the numerical solution evaluated at ``t_1``:"

# ╔═╡ 321efb40-9814-11eb-2135-dbebdbfb1e41
md"
---
### Logistic equation: $\dot{x} = f(x) := x (1 - x)$"

# ╔═╡ cf000402-9818-11eb-1924-0d651b0bf65c
@bind pow_logistic Slider(-10.0:0.20:3, default=0.6, show_value=true)

# ╔═╡ 494a5e06-9814-11eb-2d9a-095573b33364
begin
	Δt_logistic = 2.0^pow_logistic
	md"""
	#### Step size:
	
	The step size ``{\large \Delta t_{logistic} =}`` $Δt_logistic.
	
	Change it using the slider below!
	"""
end

# ╔═╡ b42a92d0-981a-11eb-04ee-d95013937cca
md"""
#### Starting from
``{\large x_0}`` = $(@bind x₀_logistic Slider(0:0.01:2, default=0.25,show_value=true))
"""

# ╔═╡ 096f350e-981e-11eb-0fe5-ed5b0ee77b8a
md"""
---
### What happens if Δt is too small?
"""

# ╔═╡ 6692d6e6-981e-11eb-2948-8f52db74a8fb
md"""
Clearly, picking a smaller Δt leads to a more accurate solution, because the truncation error at each step $n$ is smaller:

```math
x[\Delta t] - x(\Delta t) = O(\Delta t^2)
```

Does this mean we can achieve any desired degree of accuracy by taking Δt sufficiently small?

Let's investigate:
"""

# ╔═╡ f0a673ea-9825-11eb-0f18-9f5850c4f0e3
md"Here is a plot of the error (true solution - numerical solution) evaluated at the end of the period of integration ``t_1``, plotted against the step size Δt:"

# ╔═╡ 1052685c-9826-11eb-1ac6-adb8477d0e9a
md"""
Observe that as Δt decreases towards the floating point precision of the machine (around 1e-15), the error starts to increase again!

This is because floating point arithmetic introduces *round-off errors* that can be relatively significant when dealing with small numbers. As ``Δt`` decreases, the number of steps needed to span the integration interval ``(t_f - t_0)/\Delta t`` increases, and these errors accumulate.

Here are two simple examples demonstrating round-off error:
"""

# ╔═╡ 429c3e24-9820-11eb-2e54-33f9c73629f4
md"n = $(@bind n Slider(1:1000, default=49, show_value=true))"

# ╔═╡ db0bff94-981f-11eb-2f91-c7d0d5321837
begin
err1 = n * (1/n)
md"""
---
#### Does ``n \times (1/n) = 1``?


``n \times \frac{1}{n} =`` $err1
"""
end

# ╔═╡ 017ab17a-982a-11eb-0d70-758882a13532
md"""
---
#### A numerical derivative

We approximate 
```math
\frac{d}{d t} \sin(t) \approx \frac{\sin(t + \Delta t) - \sin(t)}{\Delta t}
```
and plot the deviation from the true derivative ``cos(t)``.
```math
e(Δt) := \left| \frac{\sin(t + \Delta t) - \sin(t)}{\Delta t} - \cos(t) \right|
```
"""

# ╔═╡ a45f1d9c-982b-11eb-223a-e1a07082169b
md"``{\large t = }`` $(@bind t Slider(0:0.001:2*pi,default=0.25,show_value=true))"

# ╔═╡ 654dfac0-982a-11eb-042b-41c7abaf38c9
let
	# t = 0.25
	numer_cos(Δt) = (sin(t + Δt) - sin(t))/Δt
	err2(Δt) = abs(cos(t) - numer_cos(Δt))
	plot(err2, xlims = (1e-16, 1e-8), xscale=:log10, xlabel="Δt", ylabel="e(Δt)")
end

# ╔═╡ 7d0c7dd2-982b-11eb-1bd0-c3017bd9416f
md"""
Again, the error increases drastically as ``\Delta t`` approaches the floating point precision!

Note that the above is equivalent to
```math
\sin(t + \Delta t) = \sin(t) + \cos(t) \Delta t \pm e(\Delta t) \Delta t
```
or
```math
x(t + \Delta t) = x(t) + f(x) \Delta t \pm e(\Delta t) \Delta t
```
Which is the Forward Euler integration step applied to ``\dot{x} = f(x) := \sqrt{1-x^2}``. So what we have plotted is an example of the error due to floating point calculations in the Euler method.
"""

# ╔═╡ f6cfb008-980f-11eb-04bd-c7e5aee7aea3
md"### Helper Functions"

# ╔═╡ f682b22e-9735-11eb-385d-a99c2594d12e
f_capacitor(x::Real) = 1-x

# ╔═╡ ef15b476-9813-11eb-2929-f54811a89ae1
f_logistic(x::Real) = x * (1-x)

# ╔═╡ 4bf0761a-9736-11eb-2db6-df501a9c7222
function Euler(f, x₀, (t₀,t₁), Δt)
	N = Int(floor((t₁ - t₀)/Δt))
	sol = fill(0.0, N+1, 1)
	x = x₀
	for j=1:N+1
		sol[j] = x
		x = x + f(x)*Δt
	end
	return sol
end

# ╔═╡ b075e2aa-981c-11eb-38a0-594c48470ee0
true_sol_capacitor(t) = 1 - (1 - x₀_capacitor)*exp(-t)

# ╔═╡ f72cd304-9739-11eb-16da-d95ca2fc72b8
let
	t₀, t₁ = tspan = (0.0, 5.0)
	euler_sol_capacitor = Euler(f_capacitor, x₀_capacitor, tspan, Δt_capacitor)
	# true_sol_capacitor = t -> 1 - (1 - x₀_capacitor)*exp(-t)
	
	plot(t₀:Δt_capacitor:t₁, euler_sol_capacitor, lw=2.0, alpha=0.75, lc=:red, ylims = (0.0,2.0), label="Forward Euler", xlabel="Time", ylabel="x")
	
	if (t₁-t₀ < 50 * Δt_capacitor) # only plot Scatter plot if there are a few points
		scatter!(t₀:Δt_capacitor:t₁, euler_sol_capacitor, lc=:red, label="" )
	end
	
	plot!(true_sol_capacitor, xlims = (t₀, t₁), lw=2.0, alpha=0.75, lc=:blue, label="True solution")
	
	hline!([1.0], alpha=0.6, ls=:dash, lc=:green, label="")
end

# ╔═╡ 11aeaa30-9745-11eb-1922-eb2412b05738
let
	t₀, t₁ = tspan = (0.0, 5.0)
	dtrange = range(1e-4, 1.0, length=1000)
	error = [ abs(true_sol_capacitor(t₁) - Euler(f_capacitor, x₀_capacitor, (t₀, t₁), dt)[end]) for dt in dtrange]
	plot(dtrange, error, xscale=:log10, label="Error", xlabel="Δt", ylabel="Error")
end

# ╔═╡ df47dcd4-981c-11eb-0814-9fc36dc02334
true_sol_logistic(t::Real) = x₀_logistic/(x₀_logistic + (1 - x₀_logistic)*exp(-t))

# ╔═╡ 3e247bb0-9819-11eb-03de-1d65e41e1249
let
	t₀, t₁ = tspan = (0.0, 7.0)
	
	euler_sol_logistic = Euler(f_logistic, x₀_logistic, tspan, Δt_logistic)
	
	plot(t₀:Δt_logistic:t₁, euler_sol_logistic, lw=2.0, alpha=0.75, lc=:red, ylims = (0.0,2.0), label="Forward Euler", xlabel="Time", ylabel="x")
	
	if (t₁-t₀ < 50 * Δt_logistic) # only plot Scatter plot if there are a few points
		scatter!(t₀:Δt_logistic:t₁, euler_sol_logistic, lc=:red, label="" )
	end
	
	plot!(true_sol_logistic, xlims = (t₀, t₁), lw=2.0, alpha=0.75, lc=:blue, label="True solution")
	
	hline!([0.0 1.0], ls=:dash, lc=:green, alpha=0.6, label="")
end

# ╔═╡ 6fb7519c-981d-11eb-0a03-8738d8044955
let
	t₀, t₁ = tspan = (0.0, 7.0)
	dtrange = range(1e-4, 1.0, length=1000)
	error = [ abs(true_sol_logistic(t₁) - Euler(f_logistic, x₀_logistic, (t₀, t₁), dt)[end]) for dt in dtrange]
	plot(dtrange, error, xscale=:log10, label="Error", xlabel="Δt", ylabel="Error")
end

# ╔═╡ a055b1f6-9822-11eb-38b7-63af695f266d
let
	N = 100
	t₀, t₁ = 0, 1e-12
	dtrange = 10.0.^range(-17, -14, length=100)
	error = [ abs(true_sol_logistic(t₁) - Euler(f_logistic, x₀_logistic, (t₀, t₁), dt)[end]) for dt in dtrange ]
	plot(dtrange, error, xscale=:log10, label="Error", xlabel="Δt", ylabel="Error")
end

# ╔═╡ Cell order:
# ╟─7b18828e-984e-11eb-0f5a-ff870f4b3f7e
# ╟─e43df1dc-9735-11eb-2361-8dbf881cd29d
# ╟─afac1ac4-97ef-11eb-1521-e36237bb6988
# ╟─6dd49d18-982d-11eb-3a78-8bb79e6cd97a
# ╟─75050150-97f0-11eb-10dd-45fff81be40b
# ╟─c190a43c-9739-11eb-071d-a721c779ba74
# ╟─18a72e4e-9738-11eb-30c2-a5700e9ac19f
# ╟─abdbbe84-9814-11eb-36a8-2b769e0eae65
# ╟─f72cd304-9739-11eb-16da-d95ca2fc72b8
# ╟─3c63472a-982f-11eb-361a-355745c9d948
# ╟─11aeaa30-9745-11eb-1922-eb2412b05738
# ╟─321efb40-9814-11eb-2135-dbebdbfb1e41
# ╟─494a5e06-9814-11eb-2d9a-095573b33364
# ╟─cf000402-9818-11eb-1924-0d651b0bf65c
# ╟─b42a92d0-981a-11eb-04ee-d95013937cca
# ╟─3e247bb0-9819-11eb-03de-1d65e41e1249
# ╟─6fb7519c-981d-11eb-0a03-8738d8044955
# ╟─096f350e-981e-11eb-0fe5-ed5b0ee77b8a
# ╟─6692d6e6-981e-11eb-2948-8f52db74a8fb
# ╟─f0a673ea-9825-11eb-0f18-9f5850c4f0e3
# ╟─a055b1f6-9822-11eb-38b7-63af695f266d
# ╟─1052685c-9826-11eb-1ac6-adb8477d0e9a
# ╟─db0bff94-981f-11eb-2f91-c7d0d5321837
# ╟─429c3e24-9820-11eb-2e54-33f9c73629f4
# ╟─017ab17a-982a-11eb-0d70-758882a13532
# ╟─a45f1d9c-982b-11eb-223a-e1a07082169b
# ╟─654dfac0-982a-11eb-042b-41c7abaf38c9
# ╟─7d0c7dd2-982b-11eb-1bd0-c3017bd9416f
# ╟─f6cfb008-980f-11eb-04bd-c7e5aee7aea3
# ╟─f682b22e-9735-11eb-385d-a99c2594d12e
# ╟─ef15b476-9813-11eb-2929-f54811a89ae1
# ╟─4bf0761a-9736-11eb-2db6-df501a9c7222
# ╟─b075e2aa-981c-11eb-38a0-594c48470ee0
# ╟─df47dcd4-981c-11eb-0814-9fc36dc02334
