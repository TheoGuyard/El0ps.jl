# Solving problems

In this section, we consider a [`Problem`](@ref) instance that is to be solved.

```@example optimize
using El0ps
using Random
Random.seed!(42)

y = randn(10)
f = LeastSquares(y)
M = 1.
h = Bigm(M)
A = randn(10, 30)
λ = 0.1 * compute_λmax(f, h, A)

problem = Problem(f, h, A, λ)
```

## Running the solver

The problem can be addressed with a [`BnbSolver`](@ref) that can be instantiated as follows 

```@example optimize
solver = BnbSolver()
```

Then, the [`optimize`](@ref) function allows to solve the problem.

```@example optimize
result = optimize(solver, problem)
```

The result specifies different statistics about the solving process:
* The solution status
  * `OPTIMAL` when convergence is achieved 
  * `TIME_LIMIT` when the maximum time allowed is reached
  * `ITERATION_LIMIT` when the maximum number of nodes allowed is reached
* The best objective value obtained
* The number of non-zero elements in the solution
* The last gap
* The overall solution time
* The number of nodes processed

The solution can be accessed via
```julia
result.x
```

!!! warning
    When the solution status is not `OPTIMAL`, then the objective value may not correspond to the optimal value of the problem and the solution returned may not be the optimal one.

## Specifying parameters

When creating a [`BnbSolver`](@ref), different parameters can be specified.
First, there exists parameters to control and limit the execution of the solver:
* `exploration` : the BnB exploration strategy (`BFS`, `DFS` or `MIXED`)
* `branching` : the BnB branching strategy (`LARGEST` or `RESIDUAL`)
* `maxtime` : the maximum solution time in seconds
* `maxnode` : the maximum number of nodes 
* `tolgap` : the targeted MIP duality gap in the BnB
* `tolint` : the integrality tolerance in the BnB
* `tolprune` : the tolerance when testing if a node is pruned
  
Moreover, the use can toggle different acceleration strategies using the boolean parameters. 
All these acceleration do not impact the correctness of the BnB process.
* `dualpruning` : allows to detect prunable nodes early
* `l0screening` : enable the node-screening methodology to avoid uninteresting nodes in the BnB
* `l1screening` : enable the screening methodology to accelerate the node bounding process

Finally, the BnB displays and logs can also be controlled via:
* `verbosity` : toggle displays while solving the problem
* `showevery` : difference in nodes between two consecutive displays
* `keeptrace` : return a trace of the exploration in the result

They can be passed to the solver as keywords arguments:
```@example optimize
solver = BnbSolver(maxtime=60., verbosity=false)
```

More information is given in the documentation of the [`BnbOptions`](@ref) structure that collects the keywords passed to a [`BnbSolver`](@ref).

## Warm start

When calling [`optimize`](@ref), a warm start can be specified as follows:
```@example optimize
x0 = rand(30)
optimize(solver, problem, x0=x0)
```
The BnB algorithm will construct its first upper bound based on the evaluation of `x0` in the objective of the problem.
Moreover, it is also possible to force indices of the problem variable to zero or non-zero from the beginning of the algorithm.
This can be done with
```@example optimize
force_zer = [1,2,3]
force_nnz = [4,5,6]
result = optimize(solver, problem, S0=force_zer, S1=force_nnz)
```
One notices that
```@example optimize
result.x[force_zer]
```
and
```@example optimize
result.x[force_nnz]
```
indeed correspond to zero and non-zero values, respectively.