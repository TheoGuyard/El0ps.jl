# Creating problems

This package is designed to solve L0-penalized problems.
The different components of this type of problems are: the function, $f(\cdot)$, the function $h(\cdot)$, the matrix $\mathbf{A}$ and the hyper-parameter $\lambda>0$.
Here is an example on how a [`Problem`](@ref) can be instantiated.
```@example problems
using El0ps
using Random
Random.seed!(42)  # For reproducibility

# Problem data
y = randn(10)
f = LeastSquares(y)  # Least-squares loss
M = 1.
h = Bigm(M)  # Indicator of the bound constraint -M <= x <= M
A = randn(10, 30)
λ = 0.1

# Problem instantiation
problem = Problem(f, h, A, λ);
```
The functions `f` and `h` are defined via structures provided by default in our package. 
In a [`Problem`](@ref) instance, the functions `f` and `g` must derive from an [`AbstractDatafit`](@ref) and an [`AbstractPenalty`](@ref) structure, respectively.
The number of rows in `A` must match the input dimension of `f`.
It columns needs not to be normalized but all-zero columns must be removed.
Finally, the parameter `λ` must be strictly positive.

The problem can be pretty-printed as follows:

```@example problems
println(problem)
```

The value of `λmax` is such that the all-zero vector is always solution of the problem when `λ >= λmax`.
It can be computed from the problem data as follows:
```@example problems
λmax = compute_λmax(f, h, A)
```

Note that setting `λ < λmax` does not necessarily ensure to obtain a solution with some non-zero elements.
However, we observe that it is almost always the case in practice.

!!! note 
    When the function `h` does not allows to compute `λmax` in closed form, then [`approximate_λmax`](@ref) is used to approximate `λmax` using a bisection method.

The objective value of a [`Problem`](@ref) can be evaluated as follows:
```@example problems
x = rand(30)
v = objective(problem, x)
```
If the value of `Ax` is already computed, this evaluation can also be done with
```@example problems
w = A * x
v = objective(problem, x, w)
```
in order to save computations.