# El0ps.jl

[![Build Status](https://github.com/TheoGuyard/El0ps.jl/workflows/CI/badge.svg)](https://github.com//TheoGuyard/El0ps.jl/actions)
[![Coverage](https://codecov.io/gh/TheoGuyard/El0ps.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TheoGuyard/El0ps.jl)

An Exact L0-penalized Problem Solver.

## Summary

This packages provides optimizers for the L0-penalized problem

$(\mathcal{P}) \quad \min_{\mathbf{x}} \ F(\mathbf{y},\mathbf{A}\mathbf{x}) + \lambda \|\|\mathbf{x}\|\|_0 \ \ \text{s.t.} \ \ \|\|\mathbf{x}\|\|_{\infty} \leq M$

that aim to decompose some vector $\mathbf{y} \in \mathbf{R}^{m}$ using the columns of the matrix $\mathbf{A} \in \mathbf{R}^{n \times m}$ through some model encoded into the loss function $F$.
The $\ell_0$-norm enforces sparsity in the optimizers.
The parameter $\lambda > 0$ allows to perform a trade-off between these two terms.
The Big-M constraint set with some $M > 0$ allows to devise bounded relaxations of the problem.
This problem is NP-hard and can be reformulated as a Mixed-Integer Program (MIP).

This package provides :

- A direct solution method that models $(\mathcal{P})$ as a MIP and then uses any generic solver for this type of problems (CPLEX, Gurobi, ...).
- A Branch-and-Bound method specialized for the problem $(\mathcal{P})$ enhanced with acceleration strategies leveraging the sparse structure of the problem.
- Utilities to fit a regularization path for $(\mathcal{P})$, i.e., to generate multiple solutions by varying the parameter $\lambda$.
- Utilities to generate synthetic instances of $(\mathcal{P})$.

Moreover, we provide a flexible framework to define the loss function $F$. One only has to define a function to compute its value, its gradient, its conjugate function and how to model it in a MIP.
Currently supported losses are :

- Least-squares : $F(\mathbf{y},\mathbf{w}) = \tfrac{1}{m}\|\|\mathbf{y} -\mathbf{w}\|\|_2^2$
- Logistic : $F(\mathbf{y},\mathbf{w}) = \tfrac{1}{m}\mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y} \odot \mathbf{w}))$

Please raise an issue if you want to add others.

## Installation

`El0ps.jl` can be installed using Julia's `Pkg` module.

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

To use the direct solution method, you have to install the desired MIP solvers. For instance, you have to run

```julia
pkg> add SCIP
```

in order to be able to use the [SCIP](https://github.com/scipopt/SCIP.jl) solver in the direct method.

## Examples

The directory `examples` shows classical use-cases of this package :

- `optimize.jl` : Solve a problem using the direct and the Branch-and-Bound solver.
- `path.jl` : Fit a regularization path.

## Citations

This package is linked to a paper submitted to the [ROADEF 2023](https://roadef2023.sciencesconf.org) conference.
The quote will be published soon.