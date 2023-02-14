@testset "Path" begin
    
    k = 2
    m = 30
    n = 50
    x = zeros(n)
    x[1:k] .= 1.
    A = randn(m, n)
    y = A * x
    y += 0.01 * norm(y) * randn(m) 
    f = LeastSquares(y)
    h = Bigm(1.)

    @testset "Path with BnbSolver" begin
        solver = BnbSolver(verbosity=false, maxtime=60.)
        path = fit_path(solver, f, h, A, 
            λratio_min          = 1e-2, 
            verbosity           = true,
            max_support_size    = k + 1,
            stop_if_unsolved    = true,
            compute_cv          = true,
        )
        @test all(path.converged .== true)
    end
end
