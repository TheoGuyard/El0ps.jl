# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides solution methods to address L0-penalized problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda \big(\|\mathbf{x}\|_0 + h(\mathbf{x})\big)$$

where $f(\cdot)$ is a convex and differentiable function, where $h(\cdot)$ is a convex and seprable function and where $\lambda>0$ is an hyper-parameter.

## Manual outline

```@contents
Pages = ["manual/quickstart.md", "manual/problem.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library outline

```@contents
Pages = ["library/problem.md", "library/datafit.md", "library/penalty.md", "library/solver.md", "library/path.md"]
```
