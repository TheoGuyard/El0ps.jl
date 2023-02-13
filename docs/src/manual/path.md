# Regularization path

In this section, we consider some [`Problem`](@ref) data for which we want to fit a regularization path, i.e., solving the problem for a range of `λ` values.

```@example path
using El0ps
using Random

Random.seed!(42)

f = LeastSquares()
h = Bigm(1.)
A = randn(10, 30)
y = randn(10)

solver = BnbSolver()
```

## Fitting the path

The path if fitted using the [`fit_path`](@ref) function as follows:

```@example path
path = fit_path(solver, f, h, A, y, verbosity=false)
println(path)
```

This operation returns a [`Path`](@ref) instance with several information:
* `λ/λmax`: the ratio between the current `λ` and the value `λmax` above which the solution to the problem is necessarily the all-zero vector
* `Conv`: whether the solver has converged
* `Time`: solution time
* `Fval`: value of $f(\mathbf{y},\mathbf{Ax})$
* `hval`: value of $g(\mathbf{x}) = \|\mathbf{x}\|_0 + h(\mathbf{x})$
* `Nnz`: number of non-zeros in the solution
* `CV mean`: mean cross validation error on the term $f(\mathbf{y},\mathbf{Ax})$
* `CV std`: standard deviation of the cross validation error on the term $f(\mathbf{y},\mathbf{Ax})$

## Specifying parameters

When fitting the path, the following parameters can be specified:
* `λratio_max`: the maximum value of `λ/λmax`
* `λratio_min`: the minimum value of `λ/λmax`
* `λratio_num`: the number of `λ` values in the path
* `max_support_size`: when a solution `x` is obtained in the path with `norm(x,0) > max_support_size`, the fitting is stopped
* `stop_if_unsolved`: stop the fitting when the problem at some `λ` has not been solved
* `compute_cv`: whether to compute the cross validation error
* `nb_folds`: number of folds in the cross validation evaluation
* `verbosity`: toggle displays during the fitting

They are passed to the [`fit_path`](@ref) function as keyword arguments.
```julia
path = fit_path(solver, f, h, A, y, max_support_size=5, compute_cv=false)
```

More information is given in the documentation of the [`PathOptions`](@ref) struct which handle the path parameters.