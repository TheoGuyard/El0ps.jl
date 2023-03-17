module El0ps

using Dates
using Distributions
using JuMP
using LinearAlgebra
using OrderedCollections
using Printf
using Random

version() = "v0.2"
authors() = "Theo Guyard"
contact() = "guyard.theo@gmail.com"
license() = "AGPL 3.0"
export version, authors, contact, license

include("datafit/abstract.jl")
include("datafit/leastsquares.jl")
include("datafit/logistic.jl")
include("datafit/squaredhinge.jl")
export AbstractDatafit
export dim_input, lipschitz_constant, value, gradient, conjugate, params_to_dict
export LeastSquares, Logistic, SquaredHinge

include("penalty/abstract.jl")
include("penalty/bigm.jl")
include("penalty/bigml1norm.jl")
include("penalty/bigml2norm.jl")
include("penalty/l1norm.jl")
include("penalty/l1l2norm.jl")
include("penalty/l2norm.jl")
export AbstractPenalty
export value, conjugate, prox, params_to_dict
export compute_τ, compute_μ, value_1d, conjugate_1d, prox_1d
export Bigm, BigmL1norm, BigmL2norm, L1norm, L1L2norm, L2norm

include("problem.jl")
export Problem, objective, compute_λmax

include("bounding/abstract.jl")
include("bounding/accelerations.jl")
include("bounding/cdas.jl")
export AbstractBoundingSolver
export LOWER_BOUNDING, UPPER_BOUNDING
export bound!
export CDAS

include("solver.jl")
export AbstractSolver, AbstractResult
export BFS, DFS, MIXED
export LARGEST, RESIDUAL
export OPEN, PRUNED, SOLVED, PERFECT
export BnbSolver, BnbOptions, BnbTrace, BnbResult
export optimize

include("path.jl")
export Path, PathOptions, fit_path

end
