module El0ps

using Dates
using Distributions
using JuMP
using LinearAlgebra
using Printf
using Random

version() = "v0.1"
authors() = "Theo Guyard"
contact() = "guyard.theo@gmail.com"
license() = "AGPL 3.0"

include("datafits/core.jl")
include("penalties/core.jl")
include("bounding/core.jl")
export lipschitz_constant, value, gradient, conjugate, prox, dual_scale!, params_to_dict

include("datafits/leastsquares.jl")
include("datafits/logistic.jl")
export LeastSquares, Logistic

include("penalties/bigm.jl")
include("penalties/l1norm.jl")
include("penalties/l2norm.jl")
include("penalties/l1l2norm.jl")
export Bigm, L1norm, L2norm, L1L2norm

include("problem.jl")
export Problem, objective, compute_λmax

include("solvers/core.jl")
include("solvers/bnb.jl")
include("solvers/direct.jl")
include("bounding/accelerations.jl")
include("bounding/cd.jl")
export optimize
export BnbSolver, BnbResult, BnbTrace
export BFS, DFS
export LARGEST, RESIDUAL
export DirectSolver, DirectResult

include("path.jl")
export Path, fit_path

include("data.jl")
export synthetic_data_regression, synthetic_data_classification

end