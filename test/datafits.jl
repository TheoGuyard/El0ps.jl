

@testset "Datafits" begin
    candidates = [
        (LeastSquares, "regression"),
        # (Logistic, "classification"),
    ]

    for (test_type, test_modelling) in candidates
        F = test_type()
        @testset "$F" begin
            m = 100
            w = randn(m)
            u = randn(m)
            if test_modelling == "regression"
                y = randn(m)
            elseif test_modelling == "classification"
                y = 2. .* (randn(m) .<= 0.5) .- 1.
            end
            @test lipschitz_constant(F, y) >= 0.
            @test El0ps.value(F, y, w) >= El0ps.value(F, y, u) + gradient(F, y, u)' * (w - u)
            @test El0ps.value(F, y, w) + conjugate(F, y, u) >= w' * u
        end
    end
end