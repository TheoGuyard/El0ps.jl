@testset "Problem" begin
    m = 10
    n = 30
    f = LeastSquares(randn(m))
    h = Bigm(1.)
    A = randn(m, n)
    λ = 0.1 * compute_λmax(f, h, A)
    problem = Problem(f, h, A, λ)
    println(problem)
    @test all(problem.a .≈ [norm(ai)^2 for ai in eachcol(A)])
    @test_throws AssertionError Problem(f, h, A, -1.)
    @test_throws AssertionError Problem(LeastSquares(randn(n)), h, A, λ)
end
