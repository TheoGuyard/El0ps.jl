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

This package provides solution methods to address problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda g(\mathbf{x})$$

where $g(\mathbf{x}) = \|\|\mathbf{x}\|\|_0 + h(\mathbf{x})$.
Such problems aim to minimize a loss function $f$ of some linear model $\mathbf{Ax}$.
It also enforces sparsity in the optimizers with the $\ell_0$-norm, which counts the number of non-zero entries in its input.
The function $h$ is a perturbation term required to build-up efficient numerical procedures.
In particular, this package implements a Branch-and-Bound algorithm that exploits the structure of the problem to achieve competitive performances.
It it designed to be robust to dimensionality scaling and flexible with respect to the choice of the functions $f$ and $h$.

## Project status

The package is still in beta version and is not accessible through Julia's [general registry](https://github.com/JuliaRegistries/General) yet.
An alpha version will be released soon.
The current version is tested against Julia `1.7` and `1.8` on Linux architectures.
Please report any problem on the [issue](https://github.com/TheoGuyard/El0ps.jl/issues) page.
Feel free to suggest improvements via a [pull request](https://github.com/TheoGuyard/El0ps.jl/pulls) !

## Loss and perturbation functions

Our package is designed to be flexible regarding to the expression of the functions $f$ and $h$.
The user can define his own.
It is only necessary to specify how to evaluate some of the operators associated with these functions.
The package implements by default the following functions

| Loss / Perturbation        | Expression | Parameters
|--------------|-----|---|
| Least-Squares |  $f(\mathbf{A}\mathbf{x}) = \tfrac{1}{2}\|\|\mathbf{y} - \mathbf{A}\mathbf{x}\|\|_2^2$ | Vector $\mathbf{y}$ |
| Logistic      |  $f(\mathbf{A}\mathbf{x}) = \mathbf{1}^{\top}\log(\mathbf{1} + \exp(-\mathbf{y}\odot\mathbf{A}\mathbf{x}))$ | Vector $\mathbf{y}$ |
| Big-M |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M)$ | Scalar $M > 0$ |
| Big-M + $\ell_1$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \alpha\|\|\mathbf{x}\|\|_1$ | Scalars $M,\alpha > 0$ |
| Big-M + $\ell_2$-norm      |  $h(\mathbf{x}) = \mathbb{I}(\|\|\mathbf{x}\|\|_{\infty} \leq M) + \beta\|\|\mathbf{x}\|\|_2^2$ | Scalars $M,\beta > 0$ |
| $\ell_1$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1$ | Scalar $\alpha > 0$ |
| $\ell_2$-norm      |  $h(\mathbf{x}) = \beta\|\|\mathbf{x}\|\|_2^2$ | Scalar $\beta > 0$ |
| $\ell_1\ell_2$-norm      |  $h(\mathbf{x}) = \alpha\|\|\mathbf{x}\|\|_1 + \beta\|\|\mathbf{x}\|\|_2^2$ | Scalars $\alpha,\beta > 0$ |

In the above table, $\mathbb{I}(\mathcal{C})$ denotes the convex indicator of the constraint $\mathcal{C}$ and $\odot$ denotes the Hadamard product.
Please refer to the [docs](https://theoguyard.github.io/El0ps.jl/dev) for more details.

## Citation

To cite this package, please refer to the following [paper](https://hal.science/hal-03960204/document) (in french):

```{bibtex}
@inproceedings{guyard2023solveur,
  title={Un solveur efficace pour la r{\'e}solution de probl{\`e}mes parcimonieux avec p{\'e}nalit{\'e} L0},
  author={Guyard, Theo},
  booktitle={24{\`e}me {\'e}dition du congr{\`e}s annuel de la Soci{\'e}t{\'e} Fran{\c{c}}aise de Recherche Op{\'e}rationnelle et d'Aide {\`a} la D{\'e}cision},
  year={2023}
}
```