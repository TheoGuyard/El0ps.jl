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
export version, authors, contact, license

include("datafits/core.jl")
include("penalties/core.jl")
export AbstractDatafit, AbstractPenalty
export lipschitz_constant, value, gradient, conjugate, prox, dual_scale!, params_to_dict

include("datafits/leastsquares.jl")
include("datafits/logistic.jl")
export LeastSquares, Logistic

include("penalties/zeropenalty.jl")
include("penalties/bigm.jl")
include("penalties/l1norm.jl")
include("penalties/l2norm.jl")
include("penalties/l1l2norm.jl")
include("penalties/bigml1norm.jl")
include("penalties/bigml2norm.jl")
export ZeroPenalty, Bigm, L1norm, L2norm, L1L2norm, BigmL1norm, BigmL2norm

include("problem.jl")
export Problem, objective, compute_λmax

include("solvers/core.jl")
include("bounding/core.jl")
export AbstractBoundingSolver, AbstractSolver, AbstractResult

include("solvers/bnb.jl")
include("solvers/direct.jl")
include("bounding/accelerations.jl")
include("bounding/cd.jl")
export optimize
export BnbSolver, BnbResult, BnbTrace
export BFS, DFS
export LARGEST, RESIDUAL
export DirectSolver, DirectResult
export CD

include("path.jl")
export Path, fit_path

include("data.jl")
export synthetic_data_regression, synthetic_data_classification

end
