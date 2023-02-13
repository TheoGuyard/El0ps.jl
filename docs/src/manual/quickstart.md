# Quick start

## Installation

`El0ps.jl` can be installed through Julia's [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) manager as

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

To make sure that everything went properly, you can run

```julia
pkg> test El0ps
```

The package can then be loaded with 

```@example quickstart
using El0ps
```

!!! note 
    `El0ps.jl` is tested against Julia `1.7` and `1.8`.

## Problem instantiation

The different components of problems addressed by `El0ps.jl` are the loss function $f$, the perturbation term $h$, the matrix $\mathbf{A}$, the vector $\mathbf{y}$ and the hyperparameter $\lambda$.
Here is an example on how a [`Problem`](@ref) can be instantiated.
```@example quickstart
using Random
Random.seed!(42)

# Define problem data
f = LeastSquares()
h = Bigm(1.)
A = randn(10, 30)
y = randn(10)
λ = 0.1

# Problem instantiation
problem = Problem(f, h, A, y, λ)

# Problem display
println(problem)
```
The functions `f` and `h` are defined via structures provided by our package. 
The function `f` must derive from the [`AbstractDatafit`](@ref) one and the function `h` must derive from the [`AbstractPenalty`](@ref) one.
Moreover, `A` and `y` must be a matrix and a vector with dimensions `(m,n)` and `(m,)`, respectively.
The columns of `A` needs not to be normalized but all-zero columns are not allowed.
These latter can be removed without modifying the problem solution.
Finally, the parameter `λ` must be strictly positive.

When displaying a [`Problem`](@ref), `λmax` corresponds to the value of `λ` above which the solution is always the all-zero vector.
However, setting `λ < λmax` does not necessarily ensure to obtain a solution with some non-zero elements.
The value of `λmax` can be computed from the problem data as follows:
```@example quickstart
λmax = compute_λmax(f, h, A, y)
```

The objective value of a [`Problem`](@ref) can be evaluated at any point `x` using
```@example quickstart
x = rand(30)
v = objective(problem, x)
```
If the value of `Ax` is already known, this evaluation can also be done with
```@example quickstart
w = A * x
v = objective(problem, x, w)
```
in order to avoid re-evaluating the product `Ax` and to save computations.