# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides a solution method to address problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{y},\mathbf{A}\mathbf{x}) + \lambda g(\mathbf{x})$$

where $g(x) = \|\mathbf{x}\|_0 + h(x)$.
They aims to fit a model, encoded in the loss function $f$, while forcing sparsity in the optimizers through the $\ell_0$-norm, which counts the number of non-zero entries in its argument.
The function $h$ is a perturbation term required to build-up efficient numerical procedures.
In particular, `El0ps.jl` implements a Branch-and-Bound algorithm that exploits the structure of the problem to achieve competitive performances.


## Features

* Simple problem instantiation
* Easy process to define new functions $f$ and $h$
* Branch-and-Bound algorithm with
  * Several exploration strategies
  * Several branching strategies
  * Tunable parameters
  * Efficient bounding solver
  * Structure-exploiting acceleration methods
* Routines to fit regularization paths over a range of $\lambda$

 
## Manual

```@contents
Pages = ["manual/quickstart.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library

TODO

