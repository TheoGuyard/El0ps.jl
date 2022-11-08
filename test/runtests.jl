using El0ps
using Test

using JuMP
using LinearAlgebra
using SCIP

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
include("path.jl")
include("data.jl")