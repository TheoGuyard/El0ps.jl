# El0ps.jl

*An Exact L0-penalized Problem Solver.*

[![Build Status](https://github.com/TheoGuyard/El0ps.jl/workflows/CI/badge.svg)](https://github.com//TheoGuyard/El0ps.jl/actions)
[![Coverage](https://codecov.io/gh/TheoGuyard/El0ps.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TheoGuyard/El0ps.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://theoguyard.github.io/El0ps.jl/dev)


## Summary

This packages provides routines to solve L0-penalized problems of the form

$$(\mathcal{P}) \quad \min_{\mathbf{x}} \ F(\mathbf{y},\mathbf{A}\mathbf{x}) + \lambda \|\mathbf{x}\|_0 \ \ \text{s.t.} \ \ \|\mathbf{x}\|_{\infty} \leq M$$

that aim to decompose some vector $\mathbf{y} \in \mathbf{R}^{m}$ using the columns of the matrix $\mathbf{A} \in \mathbf{R}^{n \times m}$ through some model.
The function $F$ ensures the quality of the decomposition while the $\ell_0$-norm enforces sparsity.
The parameter $\lambda > 0$ allows to perform a trade-off between these two paradigms.
The Big-M constraint set with some $M > 0$ allows to construct bounded relaxations of the problem.
This problem is NP-hard and can be reformulated as a Mixed-Integer Program (MIP).

## Solvers

This package provides :

- A **Direct** solution method that models $(\mathcal{P})$ as a MIP and then uses a generic MIP solver (CPLEX, Gurobi, ...).
- A **Branch-and-Bound** method specialized for the problem $(\mathcal{P})$ enhanced with acceleration strategies leveraging the sparse structure of the problem.
- Utilities to fit a regularization path for $(\mathcal{P})$, i.e., to generate multiple solutions by varying the parameter $\lambda$.
- Utilities to generate synthetic instances of $(\mathcal{P})$.

To use the **direct** solution method, you have to install the MIP solver you want to use. For instance, run

```julia
pkg> add SCIP
```

in order to use [SCIP](https://github.com/scipopt/SCIP.jl) to solve the MIP formulation of $(\mathcal{P})$.


## Data-fidelity functions

The data fidelity function $F$ is handled in a flexible way. Our routines only have to know how to compute its value, its gradient, its conjugate function and how to model it in a MIP.
Currently supported data-fidelity functions are :

- Least-squares : $F(\mathbf{y},\mathbf{w}) = \tfrac{1}{m}\|\mathbf{y} -\mathbf{w}\|_2^2$
- Logistic : $F(\mathbf{y},\mathbf{w}) = \tfrac{1}{m}\mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y} \odot \mathbf{w}))$

Please raise an [`issue`](https://github.com/TheoGuyard/El0ps.jl/issues) if you want to add others.
You can also make a [`pull request`](https://github.com/TheoGuyard/El0ps.jl/pulls) on your own. We recommend to mimic the definition of the [`LeastSquares`](@ref) function.
