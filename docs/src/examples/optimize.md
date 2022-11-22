## Solve a problem

> This examples shows how to construct a synthetic instance of a [`Problem`](@ref) and how to solve it using either a [`DirectSolver`](@ref) or a [`BnbSolver`](@ref).

First, import all the necessary modules.
The [SCIP](https://github.com/scipopt/SCIP.jl) optimizer will be used in the [`DirectSolver`](@ref).

```julia
using El0ps
using SCIP
```

Then, define some parameters and generate a synthetic instance.

```julia
k, m, n, ρ, s = 5, 10, 30, 0.1, 10.
x, A, y = synthetic_data_regression(k, m, n, ρ, s)
```

Instantiate a Least-squares data-fidelity function and a Big-M constraint with `M=1` as follows.

```julia
F = LeastSquares()
G = Bigm(1.)
```

Finally, compute a value `λmax` such that the all-zero vector is a solution of the [`Problem`](@ref) and set the L0-regularization strength to one tenth of this value.

```julia
λmax = compute_λmax(F, G, A, y)
λ = 0.1 * λmax
```

Now, we are ready to construct an instance of the [`Problem`](@ref).

```julia
problem = Problem(F, G, A, y, λ)
```

First, we use the **Branch-and-Bound** solver and we enable some acceleration strategies.

```julia
solver = BnbSolver(
    dualpruning = true, 
    l0screening = true, 
    l1screening = true,
)
result = optimize(solver, problem)
println(result)
```


Second, we use the **Direct** solver with the [SCIP](https://github.com/scipopt/SCIP.jl) optimizer that is called with some options.

```julia
solver = DirectSolver(
    SCIP.Optimizer, 
    options = Dict(
        "display/verblevel" => 0, 
        "limits/gap"        => 1e-4,
    )
)
result = optimize(solver, problem)
println(result)
```
