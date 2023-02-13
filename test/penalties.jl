@testset "Penalties" begin

    candidates = [
        (Zero, ()),
        (Bigm, (1.)),
        (L1norm, (1.)),
        (L2norm, (1.)),
        (L1L2norm, (1., 1.)),
        (BigmL1norm, (1., 1.)),
        (BigmL2norm, (1., 1.)),
    ]

    k, m, n, ρ, σ = 2, 10, 15, 0.1, 10.
    x, A, y = synthetic_data_classification(k, m, n, ρ, σ)
    f = LeastSquares()
    options = Dict("display/verblevel" => 0, "limits/gap" => 1e-4)
    solver = DirectSolver(SCIP.Optimizer, options=options)

    for (test_type, test_params) in candidates
        h = test_type(test_params...)
        @testset "$h utilities" begin
            x = randn(n)
            z = zeros(n)
            u = randn(m)
            r = randn(n)
            η = randn()
            λ = 0.1
            v = dual_scale!(h, A, u, λ)
            @test El0ps.value(h, x) >= 0.
            @test El0ps.value(h, z) == 0.
            @test El0ps.value(h, x) ≈ El0ps.value(h, -x)
            @test El0ps.conjugate(h, x) >= 0.
            @test El0ps.conjugate(h, z) == 0.
            @test El0ps.conjugate(h, x) ≈ El0ps.conjugate(h, -x)
            @test El0ps.value(h, x) + El0ps.conjugate(h, r) >= x' * r
            @test length(prox(h, x, η)) == length(x)
            @test all(A' * u .≈ v)
            @test El0ps.conjugate(h, v / λ) <= Inf
            @test isa(params_to_dict(h), OrderedDict)
        end

        # Do not test the bind_model!() function for penalties involving SOCP 
        # or/and SOS expressions. They require solvers that are often
        # commercial and this messes up the package CI.
        if !(typeof(h) in [Zero, L1norm, L2norm, L1L2norm, BigmL2norm])
            @testset "$h model" begin
                λmax = compute_λmax(f, h, A, y)
                problem = Problem(f, h, A, y, 0.5 * λmax)
                result = optimize(solver, problem)
                @test result.termination_status == MOI.OPTIMAL
            end
        end
    end
end
