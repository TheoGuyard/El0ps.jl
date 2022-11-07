@testset "Solver" begin
    k = 3
    m = 20
    n = 30
    ρ = 0.1
    σ = 10.
    xtrue, A, y = synthetic_data_regression(k, m, n, ρ, σ, normalize=true)
    F = LeastSquares()
    G = L2norm(1.)
    λ = 0.1 * compute_λmax(F, G, A, y)
    problem = Problem(F, G, A, y, λ)

    solver = Solver(verbosity=false, maxtime=Inf)
    result = optimize(solver, problem)
    @test result.termination_status == OPTIMAL

    tolgap = 0.1
    solver = Solver(verbosity=false, maxtime=Inf, tolgap=tolgap)
    result = optimize(solver, problem)
    @test result.relative_gap <= tolgap
end