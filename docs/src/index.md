# El0ps.jl

*An Exact L0-penalized Problem Solver.*

## Summary

This package provides solution methods to address the problem

$$\min_{\mathbf{x}} \ f(\mathbf{y},\mathbf{A}\mathbf{x}) + \lambda g(\mathbf{x})$$

where $g(x) = \|\mathbf{x}\|_0 + h(x)$.
It aims to fit $\mathbf{y}$ through a model of $\mathbf{Ax}$, encoded in the loss function $f$.
It also enforces sparsity in the optimizers with the $\ell_0$-norm, which counts the number of non-zero entries in its argument.
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
* Routines to fit regularization paths


## Citation

To cite `El0ps.jl`, please refer to the following [paper](https://hal.science/hal-03960204/document) (in french):

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
Pages = ["manual/quickstart.md", "manual/optimize.md", "manual/path.md", "manual/custom.md"]
```

## Library outline

TODO
