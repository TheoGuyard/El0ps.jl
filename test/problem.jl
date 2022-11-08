@testset "Problem" begin
    F = LeastSquares()
    G = Bigm(1.)
    x, A, y = synthetic_data_regression(5, 50, 100, 0.1, 10.)
    λ = 0.1 * compute_λmax(F, G, A, y)
    problem = Problem(F, G, A, y, λ)
    println(problem)
    @test all(problem.a .≈ [norm(ai)^2 for ai in eachcol(A)])
    @test_throws ArgumentError Problem(F, G, A, vcat(y, 1.), -1.)
    @test_throws ArgumentError Problem(F, G, A, y, -1.)
end