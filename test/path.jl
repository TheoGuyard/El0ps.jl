@testset "Path" begin

    k = 2
    m = 30
    n = 50
    x = zeros(n)
    x[1:k] .= 1.0
    A = randn(m, n)
    w = A * x
    ϵ = randn(m)
    ϵ *= norm(w, 2) / (sqrt(10.0) * norm(ϵ, 2))
    y = w + ϵ
    f = El0ps.LeastSquares(y)
    h = El0ps.Bigm(1.0)

    @testset "Path with BnbSolver" begin
        solver = El0ps.BnbSolver(verbosity = false, maxtime = 60.0)
        path = El0ps.fit_path(
            solver,
            f,
            h,
            A,
            λratio_min = 5. * 1e-2,
            verbosity = true,
            max_support_size = k + 1,
            stop_if_unsolved = true,
            compute_cv = true,
        )
        @test all(path.converged .== true)
    end
end
