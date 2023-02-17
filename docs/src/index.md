# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides solution methods to address problems of the form

$$\min_{\mathbf{x}} \ f(\mathbf{A}\mathbf{x}) + \lambda g(\mathbf{x})$$

where $g(\mathbf{x}) = \|\mathbf{x}\|_0 + h(\mathbf{x})$.
Such problems aim to minimize a loss function $f$ of some linear model $\mathbf{Ax}$.
It also enforces sparsity in the optimizers with the $\ell_0$-norm, which counts the number of non-zero entries in its input.
The function $h$ is a perturbation term required to build-up efficient numerical procedures.
In particular, this package implements a Branch-and-Bound algorithm that exploits the structure of the problem to achieve competitive performances.
It it designed to be robust to dimensionality scaling and flexible with respect to the choice of the functions $f$ and $h$.

## Features

* ⚙️ Interfaces
  * Simple problem instantiation
  * Easy process to define new functions $f$ and $h$
  * Pretty-print utilities
* 🚀 Solution methods
  * Branch-and-Bound algorithm
  * Specialized exploration strategies
  * Specialized branching strategies
  * Many tunable parameters
  * Efficient bounding solver
  * Structure-exploiting acceleration methods
  * Robust to dimension-scaling
* 📈 Utilities
  * Regularization path fitting
  * Efficiency through warm-start
  * Cross-validation computations


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

## Manual outline

```@contents
Pages = ["manual/quickstart.md", "manual/problem.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library outline

```@contents
Pages = ["library/problem.md", "library/datafits.md", "library/perturbations.md", "library/solver.md", "library/path.md"]
```
