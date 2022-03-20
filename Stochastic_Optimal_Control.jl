### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 952c2afa-4c3f-4cc4-9cfc-5c491849b250
begin 
    import Pkg 
    Pkg.activate(mktempdir()) 
    Pkg.add([
        Pkg.PackageSpec(name="JuMP"),
        Pkg.PackageSpec(name="Plots"),
        Pkg.PackageSpec(name="Ipopt"),
		Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Images"), 
		Pkg.PackageSpec(name="ImageMagick"),
		Pkg.PackageSpec(name="Colors"),
		Pkg.PackageSpec(name="PyCall"),
		Pkg.PackageSpec(name="Parameters"),  
		Pkg.PackageSpec(name="StaticArrays"),
		Pkg.PackageSpec(name="LinearAlgebra")
    ])
	using PlutoUI, Plots, JuMP, Ipopt, Images, Colors, Parameters, StaticArrays, LinearAlgebra, PyCall  
end

# ╔═╡ 3d54f452-c887-11eb-0361-27994777b10e
html"""
<div style="
position: absolute;
width: calc(100% - 30px);
border: 50vw solid #282936;
border-top: 500px solid #282936;
border-bottom: none;
box-sizing: content-box;
left: calc(-50vw + 15px);
top: -500px;
height: 500px;
pointer-events: none;
"></div>

<div style="
height: 500px;
width: 100%;
background: #282936;
color: #fff;
padding-top: 50px;
">
<span style="
font-family: Vollkorn, serif;
font-weight: 700;
font-feature-settings: 'lnum', 'pnum';
"> <p style="
font-size: 1.7rem;
opacity: .7;
"><em>UC DAVIS Intelligent Mobility Laboratory (IMobiL) </em></p>
<p style="text-align: center; font-size: 2.5rem;">
<em> Stochastic Optimal Control </em>
<p style="text-align: center; font-size: 1.5rem;">
<em> Sarder Rafee Musabbir </em>
<p style="text-align: center; font-size: 1.5rem;">
<em> Advisor: Prof. Michael Zhang </em>
</p>

<style>
body {
overflow-x: hidden;
}
</style>"""

# ╔═╡ 40a2f571-7514-4db4-b1b8-17d19b08a2df
TableOfContents(title="📚 Table of Contents", aside=true) 

# ╔═╡ b9959dc3-8e7a-489f-86e3-a7bb78f6ee21
html"<button onclick='present()'>present</button>"

# ╔═╡ a6311a60-e45e-499a-8c8c-ff0f47e682b9
md"""
1. Dynamic programming, Bellman equations, optimal value functions, value and policy iteration, shortest paths.
2. Hamilton-Jacobi-Bellman equations, approximation methods, finite horizon formulations
3. Pontryagin's maximum principle, ODE and gradient descent methods, relationship to classical mechanics.
4. Linear-quadratic-Gaussian control, Riccati equations, iterative linear approximations to nonlinear problems.
5. 


*Source: Optimal Control Theory -- Emo Todorov*

"""

# ╔═╡ 29dada5a-0425-4908-aded-0204438eb98d
md""" # Discrete Control: Bellman Equations


The **Optimal Control** problem can be formalized using **Dynamic Programming** concept as follows: find an action sequence $(u_0, u_1,...u_{n-1})$ and corresponding state sequence $(x_0, x_1,...x_{n-1})$ minimizing the total cost:

$$J(x,u) = \sum_{k=0}^{n-1} \text{cost}(x_k,u_k)$$

where, $x_{k+1} = \text{next} (x_k, u_k)$ and $u_k \in \mathcal{U}(x_k)$. The initial state $x_0 = x^{\text{init}}$ and destination state $x_n = x^{\text{dest}}$ are given. We can visualize this setting with a directed graph where the states are nodes and the actions are arrows connecting the nodes. If $\text{cost} (x, u) = 1$ for all $(x, u)$ the problem reduces to finding the shortest path from $x^{\text{init}}$ to $x^{\text{dest}}$ in the graph.

"""

# ╔═╡ 2090b1e0-bb3c-4ba7-9d5e-8e70ee2b8279
md"""
```math
\int \mathrm{e}^{-x^2} \; \mathrm{d}x
```
"""

# ╔═╡ 0e475473-4f65-4849-a3d5-28c6e38dcfec
begin
	struct Foldable{C}
		title::String
		content::C
	end

	function Base.show(io, mime::MIME"text/html", fld::Foldable)
		write(io,"<details><summary>$(fld.title)</summary><p>")
		show(io, mime, fld.content)
		write(io,"</p></details>")
	end
end

# ╔═╡ 77d23a03-e591-40ba-8b54-4f6383dc9446
Foldable("Discrete Control ?", 
md""" ## Control
The **Optimal Control** problem can be formalized using **Dynamic Programming** concept as follows: find an action sequence $(u_0, u_1,...u_{n-1})$ and corresponding state sequence $(x_0, x_1,...x_{n-1})$ minimizing the total cost:

$$J(x,u) = \sum_{k=0}^{n-1} \text{cost}(x_k,u_k)$$
	
	
""") 

# ╔═╡ c2ea5787-26d3-45d5-9406-62190ca3c9d0
begin
	
	struct TwoColumn{L, R}
		left::L
		right::R
	end

	function Base.show(io, mime::MIME"text/html", tc::TwoColumn)
		write(io, """<div style="display: flex;"><div style="flex: 50%;">""")
		show(io, mime, tc.left)
		write(io, """</div><div style="flex: 50%;">""")
		show(io, mime, tc.right)
		write(io, """</div></div>""")
	end
	
end

# ╔═╡ 804d1580-8cac-4f85-8dc0-ee1117e9893b
TwoColumn(md"Note the kink at ``x=0``!", plot(-5:5, abs))

# ╔═╡ 34262a84-0d25-4bb4-9c6b-4e3cbf80463e
@gif for φ in 0 : .1 : 2π
    plot(0 : .1 : 2π, t -> sin(t + φ))
end 

# ╔═╡ 3b81db78-3cdc-43e4-9806-106fec68ced5
let
    # show the first 10 builtin colours:
    animation = @animate for i in 1:10
        scatter([0], [0], msize=10, shape=:hexagon, mcolour=i)
    end

    gif(animation, fps=30)
end

# ╔═╡ dc9f0d30-caa5-4a0c-90cb-cac9e527acbc
@bind plot_time Slider(1:42)

# ╔═╡ d21bf93d-dbfa-480f-9b66-70b71136e303
md""" # Dynamic Programming

Dynamic Programming (DP) relies on the following obvious fact: if a given state-action sequence is optimal, and we were to remove the first state and action, the remaining sequence is also optimal (with the second state of the original sequence now acting as initial state). This is the Bellman optimality principle. 

Note the close resemblance to the Markov property of stochastic processes (a process is Markov if its future is conditionally independent of the past given the present state). 

**The optimality principle can be reworded in similar language: the choice of optimal actions in the future is independent of the past actions which led to the present state. Thus optimal state-action sequences can be constructed by starting at the final state and extending backwards. Key to this procedure is the optimal value function (or optimal cost-to-go function)**

$v(x) = \text{minimal total cost for completing the task starting from state} \quad x$


This function captures the long-term cost for starting from a given state, and makes it possible to find optimal actions through the following algorithm:

- **Consider every action available at the current state, add its immediate cost to the optimal value of the resulting next state, and choose an action for which the sum is minimal.**

- **The above algorithm is "greedy" in the sense that actions are chosen based on local information, without explicit consideration of all future scenarios. And yet the resulting actions are optimal. This is possible because the optimal value function contains all information about future scenarios that is relevant to the present choice of action.** Thus the optimal value function is an extremely useful quantity, and indeed its calculation is at the heart of many methods for optimal control.


The algorithm yields an optimal action $u = \pi (x) \in U(x)$ for every state $x$. A
mapping from states to actions is called control law or control policy. Once we have a control law $\pi : \mathcal{X} \rightarrow \mathcal{U}(\mathcal{X})$ we can start at any state $x_0$, generate action $u_0 = \pi(x_0)$, transition to state $x_1 = \text{next} (x_0, u_0)$, generate action $u_1 = \pi(x_1)$, and keep going until we reach $x^{\text{dest}}$.
Formally, an optimal control law $\pi$ satisfies,

$$\pi(x) = \arg \min_{u \in \mathcal{U}(x)} \Big\{ \text{cost} (x, u) + v (\text{next} (x, u)) \Big\}$$


The minimum in (1) may be achieved for multiple actions in the set $\mathcal{U}(x)$, which is why $\pi$ may not be unique. However the optimal value function $v$ is always uniquely defined, and satisfies,

$$v(x) = \arg \min_{u \in \mathcal{U}(x)} \Big\{ \text{cost} (x, u) + v (\text{next} (x, u)) \Big\}$$


Equations (1) and (2) are the **Bellman equations**.
If for some $x$ we already know $v (\text{next} (x, u))$ for all $u \in \mathcal{U}(x)$, then we can apply the Bellman equations directly and compute $\pi(x)$ and $v(x)$. Thus dynamic programming is particularly simple in acyclic graphs where we can start from $x^{dest}$ with $v (x^{dest}) = 0$, and perform a backward pass in which every state is visited after all its successor states have been visited. It is straightforward to extend the algorithm to the case where we are given non-zero final costs for a number of destination states (or absorbing states).


"""

# ╔═╡ f781c6ff-7414-48a3-a90b-e5791f65ff25
md"""
## Value iteration and Policy iteration
The situation is more complex in graphs with cycles. Here the Bellman equations are
still valid, but we cannot apply them in a single pass. This is because the presence of cycles makes it impossible to visit each state only after all its successors have been visited.
Instead the Bellman equations are treated as consistency conditions and used to design
iterative relaxation schemes-much like partial di⁄erential equations (PDEs) are treated as consistency conditions and solved with corresponding relaxation schemes. By "relaxation scheme" we mean guessing the solution, and iteratively improving the guess so as to make it more compatible with the consistency condition.
The two main relaxation schemes are value iteration and policy iteration. Value iteration uses only $v(x)$. We start with a guess $v^{0}$ of the optimal value function, and construct a sequence of improved guesses:

$$v^{i+1}(x) = \min_{u\in \mathcal{U}(x)} \Big\{ \text{cost}(x,u) + v^i (next (x,u)) \Big\}$$

This process is guaranteed to converge to the optimal value function v in a nite number of iterations. 

"""

# ╔═╡ 07aae1dc-16f9-4d80-a935-14421d746656
md"""
# Continuous Control: Hamilton-Jacobi-Bellman eqn.

We now turn to optimal control problems where the state $\mathbf{x} \in \mathbb{R}^{nx}$ and control $\mathbf{u} \in \mathcal{U}(x) \in \mathbb{R}^{n_u}$ are real-valued vectors. To simplify notation we will use the shortcut $\min_u$ instead of $\min_{u\in \mathcal{U}}$, although the latter is implied unless noted otherwise. consider the stochastic differential equation,

$$\text{stochastic dynamics:} \quad d\mathbf{x} = \mathbf{f}(\mathbf{x}, \mathbf{u}) dt + F (\mathbf{x}, \mathbf{u}) d\mathbf{w}$$

where, $d\mathbf{w}$ is $n_w-$dimensional Brownian motion. This is sometimes called a controlled Ito diffusion, with $f(\mathbf{x},\mathbf{u})$ being the drift and $F(\mathbf{x}, \mathbf{u})$ the diffusion coefficient. 

$\text{Meaning of the Dynamics:} \quad \mathbf{x}(t) = \mathbf{x}(0) + \int_{0}^{t} \mathbf{f}(\mathbf{x}(s), \mathbf{u}(s)) ds + \int_{0}^{t} F(\mathbf{x}(s), \mathbf{u}(s))d\mathbf{w}(s)$


$$\text{Euler Discretization:} \quad \mathbf{x}_{k+1} = \mathbf{x}_k + \Delta f(\mathbf{x}_k, \mathbf{u}_k) + \sqrt\Delta F (\mathbf{x}_k, \mathbf{u}_k)\varepsilon_k$$

where $\Delta$ is the time step, $\varepsilon_k \sim \mathcal{N}(0, I^{n_w})$ and $\mathbf{x}_k = \mathbf{x}(k \Delta)$. The $\sqrt\Delta$ term appears because
the variance of Brownian motion grows linearly with time, and thus the standard deviation of the discrete-time noise should scale as $\sqrt\Delta$.

$$\text{Cost in Continuous time:} \quad J(\mathbf{x}(\cdot), \mathbf{u}(\cdot)) = h(\mathbf{x}(t_f)) + \int_{0}^{t_f} \ell (\mathbf{x}(t), \mathbf{u}(t), t) dt$$


Objective is to find the control law $\mathbf{u} = \pi (\mathbf{x},t)$ that minimizes the expected total cost $\mathbf{J}$ for starting at a given $(\mathbf{x},t)$

$$\text{Cost in Discrete time:} \quad J(\mathbf{x.}, \mathbf{u.}) = h(\mathbf{x_n}) + \Delta \sum_{k=0}^{n-1} \ell (\mathbf{x_k}, \mathbf{u_k}, k\Delta)$$

$\text{number of time steps:} \quad n=t_f/\Delta$
"""

# ╔═╡ 37a93f54-9909-4bf7-9158-e08fac0311c6
md"""
## Derivation of the HJB equations

The state transitions are now stochastic: the probability distribution of $x_{k+1}$ given $x_k, u_k$ is the multivariate Gaussian

$$\mathbf{x}_{k+1} \sim \mathcal{N} \big( \mathbf{x}_k + \Delta f(\mathbf{x}_k, \mathbf{u}_k), \Delta S(\mathbf{x}_k, \mathbf{u}_k) \big)$$

$$S\big(\mathbf{x}, \mathbf{u} \big) = F(\mathbf{x}, \mathbf{u}) F(\mathbf{x}, \mathbf{u})^T$$

The Bellman equation for the optimal value function v is similar to $\mathbf{v(x)}$, except that $\mathbf{v}$ is now a function of space and time. 

$$v(\mathbf{x}, k) = \min_{\mathbf{u}} \Big\{ \Delta \ell(\mathbf{x}, \mathbf{u}, k\Delta) + E \big[ v (\mathbf{x} + \Delta f(\mathbf{x}, \mathbf{u}) + \xi, \ k+1) \big] \Big\}$$

$$\xi \sim \mathcal{N}(0, \Delta S(\mathbf{x}, \mathbf{u})) \quad \text{and} \quad v(\mathbf{x}, n) = h(\mathbf{x})$$

Consider the second-order Taylor-series expansion of $v$, with the time index $k+1$ suppressed for clarity:

$$v(\mathbf{x} + \delta) = v(\mathbf{x}) + \delta^T v_{\mathbf{x}}(\mathbf{x}) + \frac{1}{2} \delta^T v_{\mathbf{x}\mathbf{x}} (\mathbf{x}) \delta + o(\delta^3)$$

$$\delta = \Delta f(\mathbf{x}, \mathbf{u}) + \xi, \quad v_{\mathbf{x}}=\frac{\partial v}{\partial \mathbf{x}}, \quad v_{\mathbf{x}\mathbf{x}}=\frac{\partial^2 v}{\partial \mathbf{x} \partial \mathbf{x}}$$

Now compute the expectation of the optimal value function at the next state, using the
above Taylor-series expansion and only keeping terms up to first-order in $\Delta$. The result is:

$$E[v] = v(\mathbf{x}) + \Delta f(\mathbf{x}, \mathbf{u})^T v_{\mathbf{x}}(\mathbf{x}) + \frac{1}{2} tr \big( \Delta S(\mathbf{x}, \mathbf{u}) v_{\mathbf{x}\mathbf{x}} (\mathbf{x}) \big) + o(\Delta^2)$$


The trace term appears because,

$$E\big[ \xi^T v_{\mathbf{x}\mathbf{x}} \xi\big] = E \big[tr(\xi \xi^T v_{\mathbf{x}\mathbf{x}}) \big] = tr\big( \text{Cov} [\xi] v_{\mathbf{x}\mathbf{x}} \big) = tr(\Delta S v_{\mathbf{x}\mathbf{x}})$$

Note the **second-order** derivative $v_{\mathbf{x}\mathbf{x}}$ in the **first-order** approximation to $E[v]$. This is a recurrent theme in stochastic calculus. It is directly related to Ito's lemma, which states that if $x (t)$ is an Ito diffusion with coefficient $\sigma$, then

$$dg(x(t)) = g_x(x(t))dx(t) + \frac{1}{2} \sigma^2 g_{xx} (x(t))dt$$

Coming back to the derivation, we substitute the expression for $E [v]$ in $v(\mathbf{x}, k)$, move the term $v (x)$ outside the minimization operator (since it does not depend on u), and divide by $\Delta$. Suppressing $x, u, k$ on the right hand side, we have:

$$\frac{v(\mathbf{x},k)-v(\mathbf{x},k+1)}{\Delta} = \min_{\mathbf{u}} \Big\{ \ell + f^T v_{\mathbf{x}} + \frac{1}{2} tr(S v_{\mathbf{x}\mathbf{x}}) + o(\Delta) \Big\}$$

Recall that $t = k\Delta$, and consider the optimal value function $v (x,t)$ depned in continuous time. The left hand side in the above equation is then,

$$\frac{v(\mathbf{x},t)-v(\mathbf{x},t+\Delta)}{\Delta}$$


In the limit $\Delta \rightarrow 0$ the latter expression becomes $-\frac{\partial v}{\partial t}$, which we denote $-v_t$. Thus for $0 \leq t \leq t_f$ and $v (x, t_f) = h (x)$, the following holds:

$$v_t(\mathbf{x},t) = \min_{\mathbf{u} \in \mathcal{U}(x)} \Big\{ \ell (\mathbf{x}, \mathbf{u}, t) + f(\mathbf{x}, \mathbf{u})^T v_{\mathbf{x}} (\mathbf{x}, t) + \frac{1}{2} tr(S(\mathbf{x}, \mathbf{u}) v_{\mathbf{x}\mathbf{x}}(x,t)) \Big\}$$

Similarly to the **discrete case**, an optimal control $\pi(\mathbf{x}, t)$ is a value of $\mathbf{u}$ which achieves the minimum in the previous eqn:

$$\pi(\mathbf{x},t) = \arg \min_{\mathbf{u} \in \mathcal{U}(x)} \Big\{ \ell (\mathbf{x}, \mathbf{u}, t) + f(\mathbf{x}, \mathbf{u})^T v_{\mathbf{x}} (\mathbf{x}, t) + \frac{1}{2} tr(S(\mathbf{x}, \mathbf{u}) v_{\mathbf{x}\mathbf{x}}(x,t)) \Big\}$$

Equations $$v_t(\mathbf{x},t)$$ and $$\pi(\mathbf{x},t)$$ are the **Hamilton-Jacobi-Bellman (HJB)** equations.

"""

# ╔═╡ 2576ac29-b89c-4bc0-90ee-7eae28e5dc1b
md""" 
# Deterministic Control: Pontryagin's Maximum Principle

Optimal control theory is based on two fundamental ideas. One is dynamic programming
and the associated optimality principle, introduced by Bellman in the United States. The other is the maximum principle, introduced by Pontryagin in the Soviet Union. The maximum principle applies only to deterministic problems, and yields the same solutions as dynamic programming. Unlike dynamic programming, however, the maximum principle avoids the curse of dimensionality. Here we derive the maximum principle indirectly via the HJB equation, and directly via Lagrange multipliers. We also clarify its relationship to classical mechanics.


"""

# ╔═╡ 38422c75-9103-405c-8449-4dd76f5310c2
md"""
# Linear-Quadratic Gaussian Control

Optimal control laws can rarely be obtained in closed form. One notable exception is
the LQG case where the dynamics are linear, the costs are quadratic, and the noise (if
present) is additive Gaussian. This makes the optimal value function quadratic, and allows minimization of the Hamiltonian in closed form. Here we derive the LQG optimal controller in continuous and discrete time.

"""

# ╔═╡ 1019f72c-a305-4d0b-8d4d-5e8112b52ef3
md"""
## Continuous Case: Derivation via HJB eqn.

Consider the following stochastic optimal control problem:

$$\text{stochastic dynamics:} \quad d\mathbf{x} = (A\mathbf{x} + B\mathbf{u}) dt + F d\mathbf{w}$$

where, $d\mathbf{w}$ is $n_w-$dimensional Brownian motion. This is sometimes called a controlled Ito diffusion, with $f(\mathbf{x},\mathbf{u})$ being the drift and $F(\mathbf{x}, \mathbf{u})$ the diffusion coefficient. 


$$\text{cost rate:} \quad \ell (\mathbf{x}, \mathbf{u}) = \frac{1}{2} \mathbf{u}^{T} R \mathbf{u} + \frac{1}{2} \mathbf{x}^{T}Q\mathbf{x}$$

$$\text{final cost:} \quad h (\mathbf{x}) = \frac{1}{2} \mathbf{x}^{T}Q^f \mathbf{x}$$

where $R$ is symmetric positive-definite, $Q$ and $Q^f$ are symmetric, and $u$ is now unconstrained. We set $S = F F^{T}$ as before. The matrices $A, B, F, R, Q$ can be made time-varying without complicating the derivation below.
In order to solve for the optimal value function we will guess its parametric form, show that it satisfies the HJB equation (8), and obtain ODEs for its parameters.

$$\text{Guess Solution of the Value function:} \quad v(\mathbf{x},t) = \frac{1}{2} \mathbf{x}^T V(t) \mathbf{x} + a(t)$$

where $V (t)$ is symmetric. The boundary condition $v(x, t_f ) = h(x)$ implies $V (t_f) = Q^f$ and $a(t_f) = 0$. From $v(x,t)$ ablove we can compute the derivatives which enter into the HJB equation:

$$v_t(\mathbf{x},u) = \frac{1}{2} \mathbf{x}^T \dot V(t) \mathbf{x} + \dot a(t)$$
$$v_\mathbf{x}(\mathbf{x},t) = V(t)\mathbf{x}$$
$$v_{\mathbf{x}\mathbf{x}}(\mathbf{x},t) = V(t)$$




"""

# ╔═╡ 48277c01-21cf-4cc1-8624-cf97d411665a
md"""
Substituting these expressions in the equation below yields, 

$$-v_t(\mathbf{x},t)= \min_{u\in\mathcal{U}(\mathbf{x})} \Big\{ \ell(\mathbf{x},\mathbf{u},t) + f(\mathbf{x},\mathbf{u})^T v_\mathbf{x} (\mathbf{x},t) + \frac{1}{2} tr (S(\mathbf{x},\mathbf{u})v_{\mathbf{x}\mathbf{x}}(\mathbf{x},t)) \Big\}$$

$$-\frac{1}{2} \mathbf{x}^T \dot V(t) \mathbf{x} - \dot a(t) = \min_{\mathbf{u}} \Big\{ \frac{1}{2} \mathbf{u}^T R \mathbf{u} + \frac{1}{2} \mathbf{x}^T Q \mathbf{x} + (A\mathbf{x} + B\mathbf{u})^T V(t) \mathbf{x} + \frac{1}{2} tr(SV(t)) \Big\}$$ 

The Hamiltonian (i.e. the term inside the min operator) is quadratic in $u$, and its Hessian $R$ is positive-definite, so the optimal control can be found analytically:

$$\mathbf{u} = - R^{-1}B^T V(t) \mathbf{x}$$

With this $\mathbf{u}$, the control-dependent part of the Hamiltonian becomes,

$$\frac{1}{2} \mathbf{u}^T R \mathbf{u} + (B\mathbf{u})^T V(t) \mathbf{x} = - \frac{1}{2} \mathbf{x}^T V(t) B R^{-1} B^T V(t) \mathbf{x}$$

After grouping terms, the HJB equation reduces to,

``-\frac{1}{2} \mathbf{x}^T \dot V(t) \mathbf{x} - \dot a(t) =`` 
$$= -\frac{1}{2} \mathbf{x}^T \Big(Q + A^T V(t) + V(t) A - V(t) B R^{-1} B^T V(t)\Big) \mathbf{x} + \frac{1}{2} tr(SV(t))$$ 


where we replaced the term $2A^TV$ with $A^TV + VA$ to make the equation symmetric. This is justified because 
$\mathbf{x}^T A^TV\mathbf{x} = \mathbf{x}^T V^T A\mathbf{x} = \mathbf{x}^T VA\mathbf{x}$

Our guess of the optimal value function is correct if and only if the above equation holds for all $x$, which is the case when the $\mathbf{x}$-dependent terms are matched:

$$\text{ODE 1:} - \dot V(t) = Q + A^T V(t) + V(t) A - V(t) B R^{-1} B^TV(t)$$
$$\text{ODE 2:} -\dot a = \frac{1}{2} trace(SV(t))$$

Functions $V, a$ satisfying the above equations can obviously be found by initializing $V(t_f) = Q_f$, $a(t_f ) = 0$ and integrating the above ODEs backward in time. Thus $v(x,t)$ is the optimal value function with $V,a$ given by the ODEs, and $\mathbf{u}$ is the optimal control law (which in this case is unique).

The first line of ODEs is called a **continuous-time Riccati equation**. Note that it does not depend on the **noise covariance** $S$. Consequently the optimal control law (24) is also independent of $S$. The only effect of $S$ is on the total cost. As a corollary, the optimal control law remains the same in the deterministic case-called the linear-quadratic regulator (LQR).

"""

# ╔═╡ 58de253d-e3c6-4a48-b5cc-c9d24a371195
md"""
## Discrete Case: Derivation via Bellman eqn.

In practice one usually works with discrete-time systems. To obtain an optimal control law for the discrete-time case one could use an Euler approximation to (25), but the resulting equation is missing terms quadratic in the time step , as we will see below. Instead we apply dynamic programming directly, and obtain an exact solution to the discrete-time LQR problem. Dropping the (irrelevant) noise and discretizing the problem, we obtain

$$\text{deterministic dynamics:} \quad \mathbf{x}_{k+1} = (A\mathbf{x}_k + B\mathbf{u}_k)$$

$$\text{cost rate:} \quad \frac{1}{2} \mathbf{u}_k^{T} R \mathbf{u}_k + \frac{1}{2} \mathbf{x}_k^{T} Q \mathbf{x}_k$$

$$\text{final cost:} \quad \frac{1}{2} \mathbf{x}_n^{T}Q^f \mathbf{x}_n$$

where $n = t_f/\Delta$ and the correspondece to the continuous time problem is 

$$\mathbf{x}_k \leftarrow \mathbf{x}(k\Delta), \quad A \leftarrow (I + \Delta A), \quad B\leftarrow \Delta B,\quad R \leftarrow \Delta R, \quad Q \leftarrow \Delta Q$$

The guess for the optimal value function is again quadratic

$$v(x,k) = \frac{1}{2} \mathbf{x}^TV_k \mathbf{x}$$

with boundary condition $V_n = Q_f$. The Bellman equation (2) is

$$\frac{1}{2} \mathbf{x}^TV_k \mathbf{x} = \min_u \frac{1}{2} \Big\{ \mathbf{u}^T R \mathbf{u} + \mathbf{x}^T Q \mathbf{x} + (A\mathbf{x} + B\mathbf{u})^T V_{k+1} (A\mathbf{x} + B\mathbf{u}) \Big\}$$

As in the continuous-time case the Hamiltonian can be minimized analytically. The resulting optimal control law is,

$$\mathbf{u} = - \Big( R + B^T V_{k+1} B \Big)^{-1} B^T V_{k+1} A \mathbf{x}$$

Substituting this u in the Bellman equation, we obtain,

$$V_k = Q + A^T V_{k+1} A - A^T V_{k+1} B \Big( R + B^T V_{k+1} B \Big)^{-1} B^T V_{k+1} A$$

This completes the proof that the optimal value function is in the assumed quadratic form.

To compute $Vk$ for all $k$ we initialize $Vn = Qf$ and iterate (27) backward in time.
The optimal control law is linear in $x$, and is usually written as,

$$\mathbf{u}_k = - L_k \mathbf{x}_k$$

$$L_k = \Big( R + B^T V_{k+1} B \Big)^{-1} B^T V_{k+1} A$$

The time-varying matrix $L_k$ is called the control gain. It does not depend on the sequence of states, and therefore can be computed o› ine. Equation (27) is called a discrete-time Riccati equation. Clearly the discrete-time Riccati equation contains more terms that the continuous-time Riccati equation (25), and so the two are not identical. However one can verify that they become identical in the limit $\Delta \rightarrow 0$. To this end replace the matrices in (27) with their continuous-time analogues (26), and after rearrangement obtain,


$$\frac{V_k - V_{k+1}}{\Delta} = Q + A^TV_{k+1} + V_{k+1}A - V_{k+1}B \Big( R + \Delta B^T V_{k+1} B \Big)^{-1} B^T V_{k+1} + \frac{o(\Delta^2)}{\Delta}$$

where o absorbs terms that are second-order in $\Delta$. Taking the limit $\Delta \rightarrow 0$ yields the continuous-time Riccati equation (25).

"""

# ╔═╡ f6d7f970-e207-407c-80d4-8634bc21a37c
md"""
# Iterative LQ Gaussian (ILQG)

1. Create initial control sequence, apply it to the (nonlinear) dynamics, then obtain acorresponding state sequence.
2. Construct linear approximation to the dynamics & quadratic approximationto the cost; We get a LQG optimal control problemwith respect to the state and control deviations.
3. Solve the LQG problem, obtain optimal control deviation iteration sequence, and add it to the given control sequence. Go to step 1, or exit if converged.

"""

# ╔═╡ b80fd539-2801-48e6-b260-d57bb60bed06
md"""

1. **Create initial control sequence, obtain state sequence:**

$$\bar x_{k+1} = \bar x_k + \Delta t f(\bar x_k, \bar u_k)$$


"""

# ╔═╡ de8470ac-f2ff-4bc0-9f37-344f07d6d2fe
md"""
# Applications to Nonlinear problems
Apart from solving LQG problems, the methodology described here can be adapted to yield approximate solutions to non-LQG optimal control problems. This is done iteratively, as follows:

1. Given a control sequence, apply it to the (nonlinear) dynamics and obtain a corresponding state sequence.

2. Construct a time-varying linear approximation to the dynamics and a time-varying quadratic approximation to the cost; both approximations are centered at the statecontrol sequence obtained in step 1. This yields an LQG optimal control problem with respect to the state and control deviations.

3. Solve the resulting LQG problem, obtain the control deviation sequence, and add it to the given control sequence. Go to step 1, or exit if converged. Note that multiplying the deviation sequence by a number smaller than 1 can be used to implement linesearch.

"""

# ╔═╡ 0cc62505-c2c1-4575-a906-71a538b77b5a
md"""

Another possibility is to use differential dynamic programming (DDP), which is based
on the same idea but involves a second-order rather than a rst-order approximation to
the dynamics. In that case the approximate problem is not LQG, however one can assume
a quadratic approximation to the optimal value function and derive Riccati-like equations for its parameters. DDP and iterative LQG (iLQG) have second-order convergence in the neighborhood of an optimal solution. They can be thought of as the analog of Newton's method in the domain of optimal control. Unlike general-purpose second order methods which construct Hessian approximations using gradient information, DDP and iLQG obtain the Hessian directly by exploiting the problem structure. For deterministic problems they converge to state-control trajectories which satisfy the maximum principle, but in addition yield local feedback control laws. In our experience they are more e¢ cient that either ODE or gradient descent methods. iLQG has been generalized to stochastic systems (including multiplicative noise) and to systems subject to control constraints.

"""

# ╔═╡ 86fb3f60-e4ab-463a-b8bd-a001d96fb61d
begin  
	py"""   
	from casadi import * 
	import matplotlib.pyplot as plt

	x = MX.sym('x',2); # Two states
	p = MX.sym('p');   # Free parameter

	# Expression for ODE right-hand side
	z = 1-x[1]**2;
	rhs = vertcat(z*x[0]-x[1]+2*tanh(p),x[0])

	# ODE declaration with free parameter
	ode = {'x':x,'p':p,'ode':rhs}

	# Construct a Function that integrates over 1s
	F = integrator('F','cvodes',ode,{'tf':1})

	# Control vector
	u = MX.sym('u',4,1)

	x = [0,1]  # Initial state
	for k in range(4):
	  # Integrate 1s forward in time:
	  # call integrator symbolically
	  res = F(x0=x,p=u[k])
	  x = res["xf"]


	# NLP declaration 
	nlp = {'x':u,'f':dot(u,u),'g':x}; 

	# Solve using IPOPT
	solver = nlpsol('solver','ipopt',nlp)
	res = solver(x0=0.2,lbg=0,ubg=0)

	plt.plot(res["x"]) 
	#plt.show()
	"""
end  

# ╔═╡ b6711714-5bac-4e97-ae5e-b0a1845d584e
begin
	
	math = pyimport("math")
	math.sin(math.pi / 4)
	 
end

# ╔═╡ 040074fa-ab1d-4314-a538-8d87b02ba9e7
begin
	
	py"""
	import numpy as np

	def sinpi(x):
		return np.sin(np.pi * x)
	"""
	py"sinpi"(1)
	
end

# ╔═╡ 092c70e4-86cb-49d0-80dc-c24f94b19034
md""" # Car Race: Single Shooting
- An optimal control problem (OCP),
- solved with direct multiple-shooting.
- For more information see: [Casadi](https://web.casadi.org/)
"""

# ╔═╡ 6bfe6dca-aa48-4359-ac54-2e825c5e87a1
py"""
	from casadi import *

	N = 100 # number of control intervals

	opti = Opti() # Optimization problem

	# ---- decision variables ---------
	X = opti.variable(2,N+1) # state trajectory
	pos   = X[0,:]
	speed = X[1,:]
	U = opti.variable(1,N)   # control trajectory (throttle)
	T = opti.variable()      # final time

	# ---- objective          ---------
	opti.minimize(T) # race in minimal time

	# ---- dynamic constraints --------
	f = lambda x,u: vertcat(x[1],u-x[1]) # dx/dt = f(x,u)

	dt = T/N # length of a control interval
	for k in range(N): # loop over control intervals
	   # Runge-Kutta 4 integration
	   k1 = f(X[:,k],         U[:,k])
	   k2 = f(X[:,k]+dt/2*k1, U[:,k])
	   k3 = f(X[:,k]+dt/2*k2, U[:,k])
	   k4 = f(X[:,k]+dt*k3,   U[:,k])
	   x_next = X[:,k] + dt/6*(k1+2*k2+2*k3+k4) 
	   opti.subject_to(X[:,k+1]==x_next) # close the gaps

	# ---- path constraints -----------
	limit = lambda pos: 1-sin(2*pi*pos)/2
	opti.subject_to(speed<=limit(pos))   # track speed limit
	opti.subject_to(opti.bounded(0,U,1)) # control is limited

	# ---- boundary conditions --------
	opti.subject_to(pos[0]==0)   # start at position 0 ...
	opti.subject_to(speed[0]==0) # ... from stand-still 
	opti.subject_to(pos[-1]==1)  # finish line at position 1

	# ---- misc. constraints  ----------
	opti.subject_to(T>=0) # Time must be positive

	# ---- initial values for solver ---
	opti.set_initial(speed, 1)
	opti.set_initial(T, 1)

	# ---- solve NLP              ------
	opti.solver("ipopt") # set numerical backend
	sol = opti.solve()   # actual solve

	# ---- post-processing        ------
	from pylab import plot, step, figure, legend, show, spy

	plot(sol.value(speed),label="speed")
	plot(sol.value(pos),label="position")
	plot(limit(sol.value(pos)),'r--',label="speed limit")
	step(range(N),sol.value(U),'k',label="acceleration")
	legend(loc="upper left")

	#figure()
	#spy(sol.value(jacobian(opti.g,opti.x)))
	#figure()
	#spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))

	show()
"""

# ╔═╡ 486935c3-a65e-4f01-a21c-e42bfe805279
md"""
# Van Der Pol Oscillator: Dynamic Programming

"""

# ╔═╡ 4e644dd1-dfb5-4603-9345-fb8143b1f89d
py"""  
	from pylab import * 
	
	# End time
	T = 10.
	
	# Number of control intervals
	N = 20
	
	# Number of Runge-Kutta 4 steps per interval and step size
	NK = 20
	DT = T/(N*NK)
	
	# Number of discrete control values
	NU = 101
	
	# Number of discrete state values
	NX = 101
	
	# System dynamics, can be called with matricex
	def f(x1,x2,u):
	  x1_dot = (1 - x2*x2)*x1 - x2 + u
	  x2_dot = x1
	  q_dot  = x1*x1 + x2*x2 + u*u
	  return (x1_dot, x2_dot, q_dot)
	
	# Control enumeration
	U  = linspace(-1,1,NU)
	
	# State space enumeration
	x1 = linspace(-1,1,NX)
	x2 = linspace(-1,1,NX)
	X1,X2 = meshgrid(x1,x2)
	
	# For each control action and state, precalculate next state and stage cost
	stage_J = []
	next_x1 = []
	next_x2 = []
	for u in U:
	  # Take number of integration steps
	  X1_k = copy(X1)
	  X2_k = copy(X2)
	  Q_k = zeros(X1.shape)
	  for k in range(NK):
	    # RK4 integration for x1, x2 and q
	    k1_x1, k1_x2, k1_q = f(X1_k,                X2_k,                u)
	    k2_x1, k2_x2, k2_q = f(X1_k + DT/2 * k1_x1, X2_k + DT/2 * k1_x2, u)
	    k3_x1, k3_x2, k3_q = f(X1_k + DT/2 * k2_x1, X2_k + DT/2 * k2_x2, u)
	    k4_x1, k4_x2, k4_q = f(X1_k + DT   * k3_x1, X2_k + DT   * k3_x2, u)
	    X1_k += DT/6*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1)
	    X2_k += DT/6*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2)
	    Q_k  += DT/6*(k1_q  + 2*k2_q  + 2*k3_q  + k4_q )
	
	  # Find out which state comes next (index)
	  X1_k = matrix.round((X1_k+1)/2*(NX-1)).astype(int)
	  X2_k = matrix.round((X2_k+1)/2*(NX-1)).astype(int)
	
	  # Infinite cost if state gets out-of-bounds
	  I = X1_k  <  0; Q_k[I]=inf; X1_k[I]=0
	  I = X2_k  <  0; Q_k[I]=inf; X2_k[I]=0
	  I = X1_k >= NX; Q_k[I]=inf; X1_k[I]=0
	  I = X2_k >= NX; Q_k[I]=inf; X2_k[I]=0
	
	  # Save the stage cost and next state
	  next_x1.append(X1_k)
	  next_x2.append(X2_k)
	  stage_J.append(Q_k)
	
	# Calculate cost-to-go (no end cost) and optimal control
	J = zeros(X1.shape)
	U_opt = []
	for k in reversed(list(range(N))):
	  # Cost to go for the previous step, optimal control action
	  J_prev = inf*ones(X1.shape)
	  u_prev = -ones(X1.shape,dtype=int)
	
	  # Test all control actions
	  for uind in range(NU):
	    J_prev_test = J[next_x2[uind],next_x1[uind]]+stage_J[uind]
	    better = J_prev_test<J_prev
	    u_prev[better] = uind
	    J_prev[better] = J_prev_test[better]
	
	  # Update cost-to-go and save optimal control
	  J = J_prev
	  U_opt.append(u_prev)
	
	# Reorder U_opt by stage
	U_opt.reverse()
	
	# Find optimal control starting at x1=0, x2=1
	i1 = NX//2
	i2 = NX-1
	u_opt = []
	x1_opt = [x1[i1]]
	x2_opt = [x2[i2]]
	cost = 0
	for k in range(N):
	  # Get the optimal control and go to next step
	  u_ind = U_opt[k][i2,i1]
	  cost += stage_J[u_ind][i2,i1]
	  i1, i2 = next_x1[u_ind][i2,i1], next_x2[u_ind][i2,i1]
	
	  # Save the trajectories
	  u_opt.append(U[u_ind])
	  x1_opt.append(x1[i1])
	  x2_opt.append(x2[i2])
	
	# Optimal cost
	print("Minimal cost: ", cost)
	assert abs(cost-J[NX-1,NX//2])<1e-8 # Consistency check
	
	# Plot
	figure(1)
	clf()
	 
	# Plot optimal cost-to-go
	subplot(121)
	contourf(X1,X2,J)
	colorbar()
	xlabel('x1')
	ylabel('x2')
	title('Cost-to-go') 
	
	subplot(122)
	plot(linspace(0,T,N+1),x1_opt,'--')
	plot(linspace(0,T,N+1),x2_opt,'-.')
	step(linspace(0,T,N),u_opt,'-')
	plt.title("Dynamic programming solution")
	plt.xlabel('time')
	plt.legend(['x1 trajectory','x2 trajectory','u trajectory'])
	grid(True)
	show()
"""

# ╔═╡ 5aad7d7f-a8c2-4998-8865-4ec4201bde1d
py"""     
import numpy as np  
import control as ct  
import control.optimal as opt
import matplotlib.pyplot as plt

def vehicle_update(t, x, u, params):
    # Get the parameters for the model
    l = params.get('wheelbase', 3.)         # vehicle wheelbase
    phimax = params.get('maxsteer', 0.5)    # max steering angle (rad)

    # Saturate the steering input
    phi = np.clip(u[1], -phimax, phimax)

    # Return the derivative of the state
    return np.array([
        np.cos(x[2]) * u[0],            # xdot = cos(theta) v
        np.sin(x[2]) * u[0],            # ydot = sin(theta) v
        (u[0] / l) * np.tan(phi)        # thdot = v/l tan(phi)
    ])

def vehicle_output(t, x, u, params):
    return x                            # return x, y, theta (full state)

# Define the vehicle steering dynamics as an input/output system
vehicle = ct.NonlinearIOSystem(
    vehicle_update, vehicle_output, states=3, name='vehicle',
    inputs=('v', 'phi'), outputs=('x', 'y', 'theta'))



x0 = [0., -2., 0.]; u0 = [10., 0.]
xf = [100., 2., 0.]; uf = [10., 0.]
Tf = 10


Q = np.diag([0.1, 10, .1])    # keep lateral error low
R = np.eye(2) * 0.1
cost = opt.quadratic_cost(vehicle, Q, R, x0=xf, u0=uf)


constraints = [ opt.input_range_constraint(vehicle, [8, -0.1], [12, 0.1]) ]



horizon = np.linspace(0, Tf, 20, endpoint=True)
bend_left = [10, 0.01]        # slight left veer

result = opt.solve_ocp(
    vehicle, horizon, x0, cost, constraints, initial_guess=bend_left,
    options={'eps': 0.01})    # set step size for gradient calculation

# Extract the results
u = result.inputs
t, y = ct.input_output_response(vehicle, horizon, u, x0)
	
"""

# ╔═╡ bf84bd1c-3343-4028-b41c-9df179d11ef6
py""" 
# Plot the results
plt.subplot(3, 1, 1)
plt.plot(y[0], y[1])
plt.plot(x0[0], x0[1], 'ro', xf[0], xf[1], 'ro')
plt.xlabel("x [m]")
plt.ylabel("y [m]")

plt.subplot(3, 1, 2)
plt.plot(t, u[0])
plt.axis([0, 10, 8.5, 11.5])
plt.plot([0, 10], [9, 9], 'k--', [0, 10], [11, 11], 'k--')
plt.xlabel("t [sec]")
plt.ylabel("u1 [m/s]")

plt.subplot(3, 1, 3)
plt.plot(t, u[1])
plt.axis([0, 10, -0.15, 0.15])
plt.plot([0, 10], [-0.1, -0.1], 'k--', [0, 10], [0.1, 0.1], 'k--')
plt.xlabel("t [sec]")
plt.ylabel("u2 [rad/s]")

plt.suptitle("Lane change manuever")
plt.tight_layout()
plt.show()
"""

# ╔═╡ 9576e70e-5609-47d6-b571-9292f4016bf7
m = Model(with_optimizer(Ipopt.Optimizer, print_level=0)) 

# ╔═╡ 2e1cfa59-502b-4d25-bd6f-bda033c5a995
begin
	
	@variable(m, 0 <= p, start=1, base_name="Quantities of pizzas")
	@variable(m, 0 <= s, start=1, base_name="Quantities of sandwiches")
	@constraint(m, budget,     10p + 4s <=  80 )
	@NLobjective(m, Max, 100*p - 2*p^2 + 70*s - 2*s^2 - 3*p*s)

end

# ╔═╡ c002fa11-4a05-4137-aaf9-1a2dc307b823
begin
	optimize!(m)
	status = termination_status(m)
end

# ╔═╡ b75d340f-801d-4b70-970e-28716ee76692
with_terminal() do 
		if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(m)
		if (status == MOI.OPTIMAL)
			println("** Problem solved correctly **")
		else
			println("** Problem returned a (possibly suboptimal) solution **")
		end
		println("- Objective value : ", objective_value(m))
		println("- Optimal solutions:")
		println("pizzas: $(value.(p))")
		println("sandwitches: $(value.(s))")
		println("- Dual (budget): $(dual.(budget))")
	else
		println("The model was not solved correctly.")
		println(status)
	end
end

# ╔═╡ Cell order:
# ╟─3d54f452-c887-11eb-0361-27994777b10e
# ╠═952c2afa-4c3f-4cc4-9cfc-5c491849b250
# ╟─40a2f571-7514-4db4-b1b8-17d19b08a2df
# ╟─b9959dc3-8e7a-489f-86e3-a7bb78f6ee21
# ╟─a6311a60-e45e-499a-8c8c-ff0f47e682b9
# ╟─29dada5a-0425-4908-aded-0204438eb98d
# ╟─2090b1e0-bb3c-4ba7-9d5e-8e70ee2b8279
# ╟─7e22c38c-2a49-4b18-aeb7-645708624cf4
# ╟─f2737369-b148-4518-8eff-48cb493e4111
# ╟─0e475473-4f65-4849-a3d5-28c6e38dcfec
# ╟─77d23a03-e591-40ba-8b54-4f6383dc9446
# ╟─c2ea5787-26d3-45d5-9406-62190ca3c9d0
# ╠═804d1580-8cac-4f85-8dc0-ee1117e9893b
# ╠═34262a84-0d25-4bb4-9c6b-4e3cbf80463e
# ╠═3b81db78-3cdc-43e4-9806-106fec68ced5
# ╠═dc9f0d30-caa5-4a0c-90cb-cac9e527acbc
# ╟─d21bf93d-dbfa-480f-9b66-70b71136e303
# ╟─f781c6ff-7414-48a3-a90b-e5791f65ff25
# ╟─07aae1dc-16f9-4d80-a935-14421d746656
# ╟─37a93f54-9909-4bf7-9158-e08fac0311c6
# ╟─2576ac29-b89c-4bc0-90ee-7eae28e5dc1b
# ╟─38422c75-9103-405c-8449-4dd76f5310c2
# ╟─1019f72c-a305-4d0b-8d4d-5e8112b52ef3
# ╟─48277c01-21cf-4cc1-8624-cf97d411665a
# ╟─58de253d-e3c6-4a48-b5cc-c9d24a371195
# ╟─f6d7f970-e207-407c-80d4-8634bc21a37c
# ╟─b80fd539-2801-48e6-b260-d57bb60bed06
# ╟─de8470ac-f2ff-4bc0-9f37-344f07d6d2fe
# ╟─0cc62505-c2c1-4575-a906-71a538b77b5a
# ╠═86fb3f60-e4ab-463a-b8bd-a001d96fb61d
# ╠═b6711714-5bac-4e97-ae5e-b0a1845d584e
# ╠═040074fa-ab1d-4314-a538-8d87b02ba9e7
# ╠═092c70e4-86cb-49d0-80dc-c24f94b19034
# ╠═6bfe6dca-aa48-4359-ac54-2e825c5e87a1
# ╠═486935c3-a65e-4f01-a21c-e42bfe805279
# ╠═4e644dd1-dfb5-4603-9345-fb8143b1f89d
# ╠═5aad7d7f-a8c2-4998-8865-4ec4201bde1d
# ╠═bf84bd1c-3343-4028-b41c-9df179d11ef6
# ╠═9576e70e-5609-47d6-b571-9292f4016bf7
# ╠═2e1cfa59-502b-4d25-bd6f-bda033c5a995
# ╠═c002fa11-4a05-4137-aaf9-1a2dc307b823
# ╠═b75d340f-801d-4b70-970e-28716ee76692
