@testset "Solvers" begin
    k = 3
    m = 20
    n = 30
    ρ = 0.1
    σ = 10.
    xtrue, A, y = synthetic_data_regression(k, m, n, ρ, σ, normalize=true)
    F = LeastSquares()
    G = Bigm(1.)
    λ = 0.1 * compute_λmax(F, G, A, y)
    problem = Problem(F, G, A, y, λ)

    @testset "BnbSolver" begin
        solver = BnbSolver(verbosity=false, maxtime=60.)
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "BnbSolver with accelerations" begin
        solver = BnbSolver(
            verbosity   = false, 
            maxtime     = 60., 
            dualpruning = true,
            l0screening = true,
            l1screening = true,
        )
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "BnbSolver with verbosity" begin
        solver = BnbSolver(
            verbosity   = true, 
            maxtime     = 60., 
            dualpruning = true,
            l0screening = true,
            l1screening = true,
        )
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "DirectSolver" begin
        options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
        solver = DirectSolver(SCIP.Optimizer, options=options)
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end
end
