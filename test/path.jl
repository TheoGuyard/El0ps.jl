@testset "Path" begin
    k = 3
    m = 20
    n = 30
    ρ = 0.1
    σ = 10.
    xtrue, A, y = synthetic_data_regression(k, m, n, ρ, σ, normalize=true)
    F = LeastSquares()
    G = Bigm(1.)

    @testset "BnbSolver" begin
        solver = BnbSolver(verbosity=false, maxtime=60.)
        path = fit_path(solver, F, G, A, y, 
            λratio_min  = 1e-1, 
            verbosity   = false,
            compute_cv  = true,
        )
        @test all(path.converged .== true)
    end

    @testset "DirectSolver" begin
        options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
        solver = DirectSolver(SCIP.Optimizer, options=options)
        path = fit_path(solver, F, G, A, y, 
            λratio_min  = 1e-1, 
            verbosity   = false,
            compute_cv  = true,
        )
        @test all(path.converged .== true)
    end

    @testset "Misc" begin
        solver = BnbSolver(verbosity=false, maxtime=60.)
        path = fit_path(solver, F, G, A, y, 
            λratio_min  = 1e-1, 
            verbosity   = true,
            compute_cv  = false,
        )
        println(path)
        @test all(path.converged .== true)
    end
end