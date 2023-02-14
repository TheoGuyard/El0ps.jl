module El0ps

using Dates
using Distributions
using JuMP
using LinearAlgebra
using OrderedCollections
using Printf
using Random

version() = "v0.1"
authors() = "Theo Guyard"
contact() = "guyard.theo@gmail.com"
license() = "AGPL 3.0"
export version, authors, contact, license

include("datafits/core.jl")
include("perturbations/core.jl")
export AbstractDatafit, AbstractPerturbation
export lipschitz_constant, value, gradient, conjugate, prox, dual_scale!, params_to_dict, bind_model!

include("datafits/leastsquares.jl")
include("datafits/logistic.jl")
export LeastSquares, Logistic

include("perturbations/zero.jl")
include("perturbations/bigm.jl")
include("perturbations/l1norm.jl")
include("perturbations/l2norm.jl")
include("perturbations/l1l2norm.jl")
include("perturbations/bigml1norm.jl")
include("perturbations/bigml2norm.jl")
export Zero, Bigm, L1norm, L2norm, L1L2norm, BigmL1norm, BigmL2norm

include("problem.jl")
export Problem, objective, compute_λmax

include("solvers/core.jl")
include("bounding/core.jl")
export AbstractBoundingSolver, AbstractSolver, AbstractResult

include("solvers/bnb.jl")
include("solvers/direct.jl")
include("bounding/accelerations.jl")
include("bounding/cd.jl")
include("bounding/cdas.jl")
export optimize
export BnbSolver, BnbOptions, BnbResult, BnbTrace
export BFS, DFS
export LARGEST, RESIDUAL
export OPEN, PRUNED, SOLVED, PERFECT
export DirectSolver, DirectResult
export CD, CDAS

include("path.jl")
export Path, fit_path

include("data.jl")
export synthetic_data_regression, synthetic_data_classification

end
