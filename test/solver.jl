@testset "Solvers" begin

    k = 2
    m = 30
    n = 50
    x = zeros(n)
    x[1:k] .= 1.0
    A = randn(m, n)
    w = A * x
    y = w + (0.01 * norm(w)) * randn(m)
    f = El0ps.LeastSquares(y)
    h = El0ps.Bigm(1.0)
    λ = 0.01 * El0ps.compute_λmax(f, h, A)

    problem = El0ps.Problem(f, h, A, λ)
    maxtime = 60.0

    @testset "BnbSolver" begin
        solver = El0ps.BnbSolver(maxtime = maxtime)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "Warm start" begin
        x0 = copy(x)
        S0 = findall(x .== 0.0)
        S1 = findall(x .!= 0.0)
        solver = El0ps.BnbSolver(maxtime = maxtime)
        result = El0ps.optimize(solver, problem, x0 = x0, S0 = S0, S1 = S1)
        @test result.node_count == 1
        @test result.termination_status == MOI.OPTIMAL
        @test all(result.x[S0] .== 0.0)
        @test all(result.x[S1] .!= 0.0)
    end

    @testset "Exploration strategies" begin
        for strategy in [El0ps.BFS, El0ps.DFS, El0ps.MIXED]
            solver = El0ps.BnbSolver(maxtime = maxtime, exploration = strategy)
            result = El0ps.optimize(solver, problem)
            @test result.termination_status == MOI.OPTIMAL
        end
    end

    @testset "Branching strategies" begin
        for strategy in [El0ps.LARGEST, El0ps.RESIDUAL]
            solver = El0ps.BnbSolver(maxtime = maxtime, branching = strategy)
            result = El0ps.optimize(solver, problem)
            @test result.termination_status == MOI.OPTIMAL
        end
    end

    @testset "Node limit" begin
        solver = El0ps.BnbSolver(maxtime = maxtime, maxnode = 0)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.ITERATION_LIMIT
    end

    @testset "Time limit" begin
        solver = El0ps.BnbSolver(maxtime = 0.0)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.TIME_LIMIT
    end

    @testset "Accelerations" begin
        solver = El0ps.BnbSolver(
            maxtime = maxtime,
            dualpruning = true,
            l1screening = true,
            l0screening = true,
        )
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "Verbosity" begin
        solver = El0ps.BnbSolver(maxtime = maxtime, verbosity = true)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "Trace" begin
        solver = El0ps.BnbSolver(maxtime = maxtime, keeptrace = true)
        result = El0ps.optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end
end
