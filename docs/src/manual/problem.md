# Creating problems

This package is designed to solve L0-penalized problems.
The different components of this type of problems are: the loss function $f$, the perturbation function $h$, the feature matrix $\mathbf{A}$, and the hyperparameter $\lambda$.
Here is an example on how a [`Problem`](@ref) can be instantiated.
```@example problems
using El0ps
using Random
Random.seed!(42)

# Problem data
y = randn(10)
f = LeastSquares(y)
M = 1.
h = Bigm(M)
A = randn(10, 30)
λ = 0.1

# Problem instantiation
problem = Problem(f, h, A, λ)
```
Here, the functions `f` and `h` are defined via structures provided by default in our package. 
The function `f` must derive from the [`AbstractDatafit`](@ref) structure and the function `h` must derive from the [`AbstractPerturbation`](@ref) structure.
The number of rows in `A` must match the input dimension of `f`.
It columns needs not to be normalized but all-zero columns are not allowed.
These latter can be removed safely without modifying the problem solutions and optimal value.
Finally, the parameter `λ` must be strictly positive.

When displaying a [`Problem`](@ref), `λmax` is a value of `λ` above which the solution is always the all-zero vector.
It can be computed from the problem data as follows:
```@example problems
λmax = compute_λmax(f, h, A)
```
Note that setting `λ < λmax` does not necessarily ensure to obtain a solution with some non-zero elements.
However, we observe that it is almost always the case in practice.

The objective value of a [`Problem`](@ref) can be evaluated at any point `x` matching the dimensions of `A` using
```@example problems
x = rand(30)
v = objective(problem, x)
```
If the value of `Ax` is already known, this evaluation can also be done with
```@example problems
w = A * x
v = objective(problem, x, w)
```
in order to avoid re-evaluating the product `Ax` and to save computations.