

@testset "Datafits" begin
    candidates = [
        (LeastSquares, "regression"),
        # (Logistic, "classification"),
    ]

    for (test_type, test_modelling) in candidates
        f = test_type()
        @testset "$f" begin
            m = 100
            w = randn(m)
            u = randn(m)
            if test_modelling == "regression"
                y = randn(m)
            elseif test_modelling == "classification"
                y = 2. .* (randn(m) .<= 0.5) .- 1.
            end
            @test lipschitz_constant(f, y) >= 0.
            @test El0ps.value(f, y, w) >= El0ps.value(f, y, u) + gradient(f, y, u)' * (w - u)
            @test El0ps.value(f, y, w) + conjugate(f, y, u) >= w' * u
        end
    end
end
