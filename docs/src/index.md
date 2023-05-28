# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides solution methods to address L0-penalized problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda \big(\|\mathbf{x}\|_0 + h(\mathbf{x})\big)$$

where $f(\cdot)$ is a convex and differentiable function, where $h(\cdot)$ is a convex and seprable function and where $\lambda>0$ is an hyper-parameter.

## Quickstart

A [`Problem`](@ref) can be instantiated using either functions $f(\cdot)$ and $h(\cdot)$ already provided by the package or using user-defined ones.
It can then be solved using our [`BnbSolver`](@ref) and a regularization path can be fitted using the [`fit_path`](@ref) function.
Here is a typical workflow:
```@example optimize
using El0ps
using Random

# Data generation
m = 10
n = 30
y = rand(m)
f = LeastSquares(y)
α = 1.0
h = L2norm(α)
A = rand(m, n)
λ = 0.1 * compute_λmax(f, h, A)

# Problem instantiation
problem = Problem(f, h, A, λ)
```

```@example optimize
# Problem resolution
solver = BnbSolver(maxtime=60.)
result = optimize(solver, problem)
```

```@example optimize
# Path fitting
path = fit_path(solver, f, h, A)
```

## Manual outline

```@contents
Pages = ["manual/installation.md", "manual/problem.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library outline

```@contents
Pages = ["library/problem.md", "library/datafit.md", "library/penalty.md", "library/solver.md", "library/path.md"]
```
