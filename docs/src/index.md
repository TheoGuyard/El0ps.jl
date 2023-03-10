# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides solution methods to address L0-penalized problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda \|\mathbf{x}\|_0 + h(\mathbf{x})$$

They aim to minimize a loss function $f$ of some linear model $\mathbf{Ax}$.
It also enforces sparsity in the optimizers with the $\ell_0$-norm, which counts the number of non-zero entries in its argument.
The function $h$ allows to enforce hard or soft constraints in the problem to construct solution with better statistical properties. 
It also allows to build-up efficient numerical procedures.
In particular, this package implements a Branch-and-Bound algorithm that exploits the structure of the problem to achieve competitive performances.
It it designed to be robust to dimensionality scaling and flexible with respect to the choice of the functions $f$ and $h$.

## Manual outline

```@contents
Pages = ["manual/quickstart.md", "manual/problem.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library outline

```@contents
Pages = ["library/problem.md", "library/datafit.md", "library/penalty.md", "library/solver.md", "library/path.md"]
```
