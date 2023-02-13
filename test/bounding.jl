@testset "Bounding" begin
    k = 3
    m = 20
    n = 30
    ρ = 0.1
    σ = 10.
    xtrue, A, y = synthetic_data_regression(k, m, n, ρ, σ)
    f = LeastSquares()
    h = Bigm(1.)
    λ = 0.1 * compute_λmax(f, h, A, y)
    problem = Problem(f, h, A, y, λ)

    bounding_solvers = [
        CD(),
        CDAS(),
    ]

    for bounding_solver in bounding_solvers
        @testset "$bounding_solver without accelerations" begin
            solver = BnbSolver(
                lb_solver   = bounding_solver, 
                ub_solver   = bounding_solver, 
                verbosity   = false, 
                maxtime     = 60.,
            )
            result = optimize(solver, problem)
            @test result.termination_status == MOI.OPTIMAL
        end
        @testset "$bounding_solver with accelerations" begin
            solver = BnbSolver(
                lb_solver   = bounding_solver, 
                ub_solver   = bounding_solver, 
                dualpruning = true,
                l0screening = true,
                l1screening = true,
                verbosity   = false, 
                maxtime     = 60.,
            )
            result = optimize(solver, problem)
            @test result.termination_status == MOI.OPTIMAL
        end
    end
end
