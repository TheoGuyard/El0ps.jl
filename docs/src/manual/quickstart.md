# Quick start

## Installation

`El0ps.jl` can be installed through Julia's [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) as follows:

```julia
pkg> add "https://github.com/TheoGuyard/El0ps.jl"
```

Then, the package is loaded with 

```@example quickstart
using El0ps
```

To make sure that it is install properly, you can run the following command:

```julia
pkg> test El0ps
```

## Problem instantiation

The ingredients of problems addressed by `El0ps.jl` are the loss function $f$, the perturbation term $h$, the matrix $\mathbf{A}$, the vector $\mathbf{y}$ and the hyperparameter $\lambda$.
Here is an example on how a problem can be instantiated.
```@example quickstart
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
In the display, `λmax` corresponds to the value of `λ` above which the solution of the problem is the all-zero vector.
Displaying the problem tells how close is the chosen `λ` from `λmax`.
The value of `λmax` can be computed from the problem data as 
```@example quickstart
λmax = compute_λmax(f, h, A, y)
```
