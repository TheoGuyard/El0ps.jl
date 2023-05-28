# El0ps.jl

*An Exact L0-penalized Problem Solver.*

[![Build Status](https://github.com/TheoGuyard/El0ps.jl/workflows/CI/badge.svg)](https://github.com//TheoGuyard/El0ps.jl/actions)
[![Coverage](https://codecov.io/gh/TheoGuyard/El0ps.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TheoGuyard/El0ps.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://theoguyard.github.io/El0ps.jl/dev)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Summary

This package provides solution methods to address L0-penalized problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda \big(\|\|\mathbf{x}\|\|_0 + h(\mathbf{x})\big)$$

where $f(\cdot)$ is a convex and differentiable function, where $h(\cdot)$ is a convex and seprable function and where $\lambda>0$ is an hyper-parameter.


## Features

**🖥 Interfaces**
* Simple problem instantiation
* Allows user-defined functions $f(\cdot)$ and $h(\cdot)$
* Popular choices of $f(\cdot)$ and $h(\cdot)$ already implemented

* Regularization path fitting
  
**🚀 Solution methods**
* Branch-and-Bound algorithm
* Specialized exploration strategies
* Specialized branching strategies
* Efficient bounding solver
* Structure-exploiting acceleration methods
* Robust to dimensionality scaling


## Installation

Our package can be installed through Julia's `Pkg` manager as follows:

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

## Project status

The package is still in beta version and is not accessible through Julia's [general registry](https://github.com/JuliaRegistries/General) yet.
An alpha version will be released soon.
The current version is tested against Julia `1.7` and `1.8` on Linux architectures.
Please report any problem on the [issue](https://github.com/TheoGuyard/El0ps.jl/issues) page.
Feel free to suggest improvements via a [pull request](https://github.com/TheoGuyard/El0ps.jl/pulls) !

## License

This software is distributed under the [GNU AGPL-3](https://www.gnu.org/licenses/agpl-3.0.en.html) license.

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
