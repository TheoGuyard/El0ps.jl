# El0ps.jl

*An Exact L0-penalized Problem Solver.*

[![Build Status](https://github.com/TheoGuyard/El0ps.jl/workflows/CI/badge.svg)](https://github.com//TheoGuyard/El0ps.jl/actions)
[![Coverage](https://codecov.io/gh/TheoGuyard/El0ps.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TheoGuyard/El0ps.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://theoguyard.github.io/El0ps.jl/dev)

## Installation

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

## Summary

This package provides a solution method to address problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda g(\mathbf{x})$$

where $g(x) = \|\|\mathbf{x}\|\|_0 + h(x)$.
They aims to fit a model, encoded in the loss function $f$, while forcing sparsity in the optimizers through the $\ell_0$-norm, which counts the number of non-zero entries in its argument.
The function $h$ is a perturbation term required to build-up efficient numerical procedures.
In particular, `El0ps.jl` implements a Branch-and-Bound algorithm that exploits the structure of the problem to achieve competitive performances.


## Loss and perturbation functions

Our implementation is designed to be flexible regarding to the expression of the functions $f$ and $h$.
The user can define his own.
It is only necessary to indicate how to evaluate some of the operators associated with them (see the [docs](https://theoguyard.github.io/El0ps.jl/dev) for more details).
The package already supports the following function $f$ and $g$.

| Loss / Perturbation        | Expression | Parameters
|--------------|-----|---|
| Least-Squares |  $f(\mathbf{A}\mathbf{x}) = \tfrac{1}{2}\|\|\mathbf{y} - \mathbf{A}\mathbf{x}\|\|_2^2$ | None |
| Logistic      |  $f(\mathbf{A}\mathbf{x}) = \mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y}\odot\mathbf{A}\mathbf{x}))$ | None |
| Big-M |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M)$ | $M > 0$ |
| Big-M + $\ell_1$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \alpha\|\|\mathbf{x}\|\|_1$ | $M,\alpha > 0$ |
| Big-M + $\ell_2$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \beta\|\|\mathbf{x}\|\|_2^2$ | $M,\beta > 0$ |
| $\ell_1$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1$ | $\alpha > 0$ |
| $\ell_2$-norm      |  $h(\mathbf{x}) = \beta\|\|\mathbf{x}\|\|_2^2$ | $\beta > 0$ |
| $\ell_1\ell_2$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1 + \beta\|\|\mathbf{x}\|\|_2^2$ | $\alpha,\beta > 0$ |

In the above table, the function $\mathbb{I}(\mathcal{C})$ denotes the convex indicator of the constraint $\mathcal{C}$.
Please raise an [issue](https://github.com/TheoGuyard/El0ps.jl/issues) or create a [pull request](https://github.com/TheoGuyard/El0ps.jl/pulls) to ask for the support of other loss and perturbation functions.

## Citation

To cite `El0ps.jl`, please refer to the following [paper](https://hal.science/hal-03960204/document) (in french) :

```{bibtex}
@inproceedings{guyard2023solveur,
  title={Un solveur efficace pour la r{\'e}solution de probl{\`e}mes parcimonieux avec p{\'e}nalit{\'e} L0},
  author={Guyard, Theo},
  booktitle={24{\`e}me {\'e}dition du congr{\`e}s annuel de la Soci{\'e}t{\'e} Fran{\c{c}}aise de Recherche Op{\'e}rationnelle et d'Aide {\`a} la D{\'e}cision},
  year={2023}
}
```