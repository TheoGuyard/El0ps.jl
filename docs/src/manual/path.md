# Regularization path

In this section, we consider a [`Problem`](@ref) for which we want to fit a regularization path, i.e., solve it over a range of `λ` values.
For each value, the corresponding problem will be solved with a [`BnbSolver`](@ref).

```@example path
using El0ps
using Random
Random.seed!(42)

y = randn(10)
f = LeastSquares(y)
M = 1.
h = Bigm(M)
A = randn(10, 30)

solver = BnbSolver()
```

## Fitting the path

The path if fitted using the [`fit_path`](@ref) function as follows:

```@example path
path = fit_path(solver, f, h, A)
```

This operation returns a [`Path`](@ref) instance with several information:
* `λ`: the current value of `λ` in the path
* `λ/λmax`: the value of the ratio `λ/λmax`
* `x`: the solution of the problem at the current value of `λ`
* `converged`: whether the solver has converged
* `solve_time`: solve time
* `node_count`: number of nodes explored
* `objective_value`: optimal objective value
* `datafit_value`: value of $f(\cdot)$ at optimum
* `penalty_value`: value of $h(\cdot)$ at optimum
* `support_size`: number of non-zeros in the solution
* `cv_mean`: mean cross validation error on the term $f(\cdot)$
* `cv_std`: standard deviation of the cross validation error on the term $f(\cdot)$

## Specifying parameters

When fitting the path, the following parameters can be specified:
* `λratio_max`: the maximum value of `λ/λmax`
* `λratio_min`: the minimum value of `λ/λmax`
* `λratio_num`: the number of `λ` values in the path
* `max_support_size`: stop the fitting when a solution `x` with `norm(x,0) > max_support_size` is obtained
* `stop_if_unsolved`: stop the fitting when the problem at some `λ` has not been solved
* `compute_cv`: whether to compute the cross validation error
* `nb_folds`: number of folds in the cross validation evaluation
* `verbosity`: toggle displays during the fitting

They are passed to the [`fit_path`](@ref) function as keyword arguments.
```julia
path = fit_path(solver, f, h, A, max_support_size=5, compute_cv=false)
```

More information is given in the documentation of the [`PathOptions`](@ref) struct which handle the path parameters.