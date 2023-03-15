@testset "Datafit" begin

    @testset "Abstract" begin
        struct NewF <: El0ps.AbstractDatafit end
        f = NewF()
        x = randn(10)
        @test_throws ErrorException El0ps.dim_input(f)
        @test_throws ErrorException El0ps.lipschitz_constant(f)
        @test_throws ErrorException El0ps.value(f, x)
        @test_throws ErrorException El0ps.gradient(f, x)
        @test_throws ErrorException El0ps.conjugate(f, x)
    end

    @testset "Candidates" begin
        candidates = [
            (El0ps.LeastSquares, m -> randn(m)),
            (El0ps.Logistic, m -> 2.0 .* (randn(m) .<= 0.5) .- 1.0),
        ]
        for (test_type, y_generation) in candidates
            @testset "$test_type" begin
                m = 100
                y = y_generation(m)
                w = randn(m)
                u = randn(m)
                f = test_type(y)
                @test isa(println(f), Nothing)
                @test El0ps.dim_input(f) == m
                @test El0ps.lipschitz_constant(f) >= 0.0
                @test El0ps.value(f, w) >= El0ps.value(f, u) + El0ps.gradient(f, u)' * (w - u)
                if !isa(f, El0ps.Logistic)
                    @test El0ps.value(f, w) + El0ps.conjugate(f, u) >= w' * u
                end
                @test isa(El0ps.params_to_dict(f), OrderedDict)
            end
        end
    end
end
