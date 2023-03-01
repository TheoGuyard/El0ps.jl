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

include("datafits/abstract.jl")
include("datafits/leastsquares.jl")
include("datafits/logistic.jl")
export AbstractDatafit
export dim_input, lipschitz_constant, value, gradient, conjugate, params_to_dict
export LeastSquares, Logistic

include("perturbations/abstract.jl")
include("perturbations/bigm.jl")
include("perturbations/bigml1norm.jl")
include("perturbations/bigml2norm.jl")
include("perturbations/l1norm.jl")
include("perturbations/l1l2norm.jl")
include("perturbations/l2norm.jl")
export AbstractPerturbation
export compute_τ,
    compute_μ, value_1d, value, conjugate_1d, conjugate, prox_1d, prox, params_to_dict
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
