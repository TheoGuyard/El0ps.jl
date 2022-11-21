@testset "Penalties" begin

    candidates = [
        (ZeroPenalty, ()),
        (Bigm, (1.)),
        (L1norm, (1.)),
        (L2norm, (1.)),
        (L1L2norm, (1., 1.)),
        (BigmL1norm, (1., 1.)),
        (BigmL2norm, (1., 1.)),
    ]

    k, m, n, ρ, σ = 2, 10, 15, 0.1, 10.
    x, A, y = synthetic_data_classification(k, m, n, ρ, σ, normalize=true)
    F = LeastSquares()
    options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
    solver = DirectSolver(SCIP.Optimizer, options=options)

    for (test_type, test_params) in candidates
        G = test_type(test_params...)
        @testset "$G utilities" begin
            x = randn(n)
            z = zeros(n)
            v = randn(n)
            η = randn()
            @test El0ps.value(G, x) >= 0.
            @test El0ps.value(G, z) == 0.
            @test El0ps.value(G, x) ≈ El0ps.value(G, -x)
            @test conjugate(G, x) >= 0.
            @test conjugate(G, z) == 0.
            @test conjugate(G, x) ≈ conjugate(G, -x)
            @test El0ps.value(G, x) + conjugate(G, v) >= x' * v
            @test length(prox(G, x, η)) == length(x)
            @test isa(params_to_dict(G), Dict)
        end

        # Do not test the bind_model!() function for penalties involving SOCP 
        # or/and SOS expressions. They require solvers that are often
        # commercial and this messes up the package CI.
        if !(typeof(G) in [ZeroPenalty, L1norm, L2norm, L1L2norm, BigmL2norm])
            @testset "$G model" begin
                λmax = compute_λmax(F, G, A, y)
                problem = Problem(F, G, A, y, 0.5 * λmax)
                result = optimize(solver, problem)
                @test result.termination_status == MOI.OPTIMAL
            end
        end
    end
end
