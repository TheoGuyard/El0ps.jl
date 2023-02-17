@testset "Problem" begin
    m = 10
    n = 30
    f = El0ps.LeastSquares(randn(m))
    h = El0ps.Bigm(1.)
    A = randn(m, n)
    λ = 0.1 * El0ps.compute_λmax(f, h, A)
    problem = El0ps.Problem(f, h, A, λ)
    @test isa(println(problem), Nothing)
    @test all(problem.a .≈ [norm(ai)^2 for ai in eachcol(A)])
    @test_throws AssertionError El0ps.Problem(f, h, A, -1.)
    @test_throws AssertionError El0ps.Problem(El0ps.LeastSquares(randn(n)), h, A, λ)
end
