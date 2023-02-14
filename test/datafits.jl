@testset "Datafits" begin
    candidates = [
        (LeastSquares, m -> randn(m)),
        (Logistic, m -> 2. .* (randn(m) .<= 0.5) .- 1.),
    ]
    for (test_type, y_generation) in candidates
        @testset "$test_type" begin
            m = 100
            y = y_generation(m)
            w = randn(m)
            u = randn(m)
            f = test_type(y)
            @test El0ps.dim_input(f) == m
            @test El0ps.lipschitz_constant(f) >= 0.
            @test El0ps.value(f, w) >= El0ps.value(f, u) + El0ps.gradient(f, u)' * (w - u)
            if !isa(f, Logistic)
                @test El0ps.value(f, w) + El0ps.conjugate(f, u) >= w' * u
            end
            @test isa(El0ps.params_to_dict(f), OrderedDict)
        end
    end
end
