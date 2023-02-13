@testset "Solvers" begin
    k = 3
    m = 20
    n = 30
    ρ = 0.1
    σ = 10.
    x, A, y = synthetic_data_regression(k, m, n, ρ, σ)
    f = LeastSquares()
    h = Bigm(1.)
    λ = 0.1 * compute_λmax(f, h, A, y)
    problem = Problem(f, h, A, y, λ)

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

    @testset "BnbSolver with trace" begin
        solver = BnbSolver(
            verbosity   = false, 
            maxtime     = 60., 
            dualpruning = true,
            l0screening = true,
            l1screening = true,
            keeptrace   = true,
        )
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "BnbSolver with initial S0 and S1" begin
        solver = BnbSolver(
            verbosity   = false, 
            maxtime     = 60., 
            dualpruning = true,
            l0screening = true,
            l1screening = true,
        )
        S0 = findall(x .== 0.)
        S1 = findall(x .!= 0.)
        result = optimize(solver, problem, S0=S0, S1=S1)
        @test result.node_count == 1
        @test result.termination_status == MOI.OPTIMAL
        @test all(result.x[S0] .== 0.)
        @test all(result.x[S1] .!= 0.)
    end

    @testset "DirectSolver" begin
        options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
        solver = DirectSolver(SCIP.Optimizer, options=options)
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "DirectSolver with initial S0 and S1" begin
        options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
        solver = DirectSolver(SCIP.Optimizer, options=options)
        S0 = findall(x .== 0.)
        S1 = findall(x .!= 0.)
        result = optimize(solver, problem, S0=S0, S1=S1)
        @test result.node_count == 1
        @test result.termination_status == MOI.OPTIMAL
        @test all(result.x[S0] .== 0.)
        @test all(result.x[S1] .!= 0.)
    end
end
