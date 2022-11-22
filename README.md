# El0ps.jl

*An Exact L0-penalized Problem Solver.*

[![Build Status](https://github.com/TheoGuyard/El0ps.jl/workflows/CI/badge.svg)](https://github.com//TheoGuyard/El0ps.jl/actions)
[![Coverage](https://codecov.io/gh/TheoGuyard/El0ps.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TheoGuyard/El0ps.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://theoguyard.github.io/El0ps.jl/dev)

## Summary

This packages provides routines to solve L0-penalized problems of the form

$$(\mathcal{P}) \quad \min_{\mathbf{x}} \ F(\mathbf{y},\mathbf{A}\mathbf{x}) + \lambda \|\|\mathbf{x}\|\|\_0 \ \ \text{s.t.} \ \ \|\|\mathbf{x}\|\|\_{\infty} \leq M$$

that aim to decompose some vector $\mathbf{y} \in \mathbf{R}^{m}$ using the columns of the matrix $\mathbf{A} \in \mathbf{R}^{n \times m}$ through some model.
The function $F$ ensures the quality of the decomposition while the $\ell_0$-norm enforces sparsity.
The parameter $\lambda > 0$ allows to perform a trade-off between these two paradigms.
The Big-M constraint set with some $M > 0$ allows to construct bounded relaxations of the problem.
This problem is NP-hard and can be reformulated as a Mixed-Integer Program (MIP).

## Installation

This package is tested against Julia `1.7` and `1.8`. It can be installed as follows.

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

## Examples

Use-cases are provided in the [documentation](https://theoguyard.github.io/El0ps.jl/dev).

## Citations

This package is linked to a paper submitted to the 2023 [ROADEF](https://www.roadef.org/) conference.
The citation will be added soon.
