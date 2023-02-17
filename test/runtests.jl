using El0ps
using Test

using JuMP
using LinearAlgebra
using OrderedCollections
using Statistics

@testset "Package info" begin
    println(version())
    println(authors())
    println(contact())
    println(license())
    @test true
end

include("datafits.jl")
include("perturbations.jl")
include("problem.jl")
include("solver.jl")
include("path.jl")
