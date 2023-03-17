@testset "Problem" begin

    @testset "Problem utilities" begin
        m = 10
        n = 30
        f = El0ps.LeastSquares(randn(m))
        h = El0ps.Bigm(1.0)
        A = randn(m, n)
        v = norm(A' * gradient(f, zeros(dim_input(f))), Inf)
        λmax = El0ps.compute_λmax(f, h, A)
        problem = El0ps.Problem(f, h, A, 0.1 * λmax)
        @test isa(println(problem), Nothing)
        @test all(problem.a .≈ [norm(ai)^2 for ai in eachcol(A)])
        @test_throws AssertionError El0ps.Problem(f, h, A, -1.0)
        @test_throws AssertionError El0ps.Problem(
            El0ps.LeastSquares(randn(n)),
            h,
            A,
            0.1 * λmax,
        )
    end

    @testset "λmax" begin
        m = 10
        n = 30
        f = El0ps.LeastSquares(randn(m))
        h = El0ps.Bigm(1.0)
        A = randn(m, n)
        v = norm(A' * gradient(f, zeros(dim_input(f))), Inf)
        λmax = El0ps.compute_λmax(f, h, A)
        problem = El0ps.Problem(f, h, A, λmax)
        solver = El0ps.BnbSolver(maxtime = 60.0)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
        @test all(result.x .≈ 0.0)
    end
end
