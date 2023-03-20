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

@testset "El0ps" begin
    include("datafit.jl")
    include("penalty.jl")
    include("problem.jl")
    include("bounding.jl")
    include("solver.jl")
    include("path.jl")
    include("utils.jl")
end
