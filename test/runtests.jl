using El0ps
using Test

using JuMP
using LinearAlgebra
using OrderedCollections
using SCIP
using Statistics

@testset "Package info" begin
    println(version())
    println(authors())
    println(contact())
    println(license())
    @test true
end

include("datafits.jl")
include("penalties.jl")
include("problem.jl")
include("solvers.jl")
include("bounding.jl")
include("path.jl")
include("data.jl")
