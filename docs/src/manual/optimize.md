# Solving problems

In this section, we consider a [`Problem`](@ref) instance that is to be solved.

```@example optimize
using El0ps
using Random

Random.seed!(42)

f = LeastSquares()
h = Bigm(1.)
A = randn(10, 30)
y = randn(10)
λ = 0.1

problem = Problem(f, h, A, y, λ)
```

## Running the solver

The problem can be addressed with the [`BnbSolver`](@ref), instantiated with 

```@example optimize
solver = BnbSolver()
```

Here, we do not specify any argument, so the default ones are used.
Then, simply call the [`optimize`](@ref) function to solve the problem.

```@example optimize
result = optimize(solver, problem)
```

The result specifies different statistics about the solution process:
* The solution status
  * `OPTIMAL` when convergence is met 
  * `TIME_LIMIT` when the maximum time allowed is reached
  * `ITERATION_LIMIT` when the maximum number of nodes allowed is reached
* The best objective value obtained
* The number of non-zero elements in the solution of the problem
* The last MIP gap, which should be zero at optimality
* The overall solution time
* The number of nodes processed

The solution of the problem can be accessed via
```@example optimize
result.x
```

!!! warning
    When the solution status is not `OPTIMAL`, then the objective value may not correspond to the optimal value of the problem and the solution returned may not be the optimal one.

## Specifying parameters

When creating a [`BnbSolver`](@ref), different parameters can be specified.
First, there exists parameters to control and limit the execution of the solver:
* `exploration` : the BnB exploration strategy (`BFS` or `DFS`)
* `branching` : the BnB branching strategy (`LARGEST` or `RESIDUAL`)
* `maxtime` : the maximum solution time allowed
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

More information is given in the documentation of the [`BnbOptions`](@ref) structure that collects the keywords passed to [`BnbSolver`](@ref).

## Warm start

When calling [`optimize`](@ref), a warm start `x0` can be specified as follows:
```@example optimize
x0 = randn(30)
result = optimize(solver, problem, x0=x0)
```
The BnB algorithm will construct its first upper bound based on the evaluation of `x0` in the objective of the problem.