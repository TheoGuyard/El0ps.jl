@testset "Solvers" begin
    
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
    λ = 0.1 * compute_λmax(f, h, A)
    problem = Problem(f, h, A, λ)

    @testset "BnbSolver" begin
        solver = BnbSolver(maxtime=60., verbosity=false)
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "BnbSolver with arguments" begin
        solver = BnbSolver(
            maxtime     = 60.,
            verbosity   = true, 
            dualpruning = true,
            l0screening = true,
            l1screening = true,
            keeptrace   = true,
        )
        result = optimize(solver, problem)
        @test result.termination_status == MOI.OPTIMAL
    end

    @testset "BnbSolver with warm start" begin
        solver = BnbSolver(
            maxtime     = 60.,
            verbosity   = false,
        )
        x0 = copy(x)
        S0 = findall(x .== 0.)
        S1 = findall(x .!= 0.)
        result = optimize(solver, problem, x0=x0, S0=S0, S1=S1)
        @test result.node_count == 1
        @test result.termination_status == MOI.OPTIMAL
        @test all(result.x[S0] .== 0.)
        @test all(result.x[S1] .!= 0.)
    end
end
