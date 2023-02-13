@testset "Problem" begin
    f = LeastSquares()
    h = Bigm(1.)
    x, A, y = synthetic_data_regression(5, 50, 100, 0.1, 10.)
    λ = 0.1 * compute_λmax(f, h, A, y)
    problem = Problem(f, h, A, y, λ)
    println(problem)
    @test all(problem.a .≈ [norm(ai)^2 for ai in eachcol(A)])
    @test_throws AssertionError Problem(f, h, A, vcat(y, 1.), -1.)
    @test_throws AssertionError Problem(f, h, A, y, -1.)
end
