@testset "Penalty" begin

    @testset "Abstract" begin
        struct NewH <: El0ps.AbstractPenalty end
        f = El0ps.LeastSquares(randn(1))
        h = NewH()
        A = randn(1, 1)
        x = randn()
        η = rand()
        λ = rand()
        @test_throws ErrorException El0ps.compute_τ(h, λ)
        @test_throws ErrorException El0ps.compute_μ(h, λ)
        @test_throws ErrorException El0ps.compute_λmax(f, h, A)
        @test_throws ErrorException El0ps.value_1d(h, x)
        @test_throws ErrorException El0ps.conjugate_1d(h, x)
        @test_throws ErrorException El0ps.prox_1d(h, x, η)
    end

    @testset "Candidates" begin
        candidates = [
            (El0ps.Bigm, (1.0)),
            (El0ps.L1norm, (1.0)),
            (El0ps.L2norm, (1.0)),
            (El0ps.L1L2norm, (1.0, 1.0)),
            (El0ps.BigmL1norm, (1.0, 1.0)),
            (El0ps.BigmL2norm, (1.0, 1.0)),
        ]
        for (test_type, test_params) in candidates
            @testset "$test_type utilities" begin
                m = 50
                n = 100
                f = El0ps.LeastSquares(randn(m))
                A = randn(m, n)
                h = test_type(test_params...)
                λ = rand()
                τ = El0ps.compute_τ(h, λ)
                μ = El0ps.compute_μ(h, λ)
                x = randn(n)
                z = zeros(n)
                r = randn(n)
                η = randn()
                λmax = El0ps.compute_λmax(f, h, A)
                @test isa(println(h), Nothing)
                if μ < Inf
                    @test El0ps.conjugate_1d(h, τ) ≈ λ
                    @test El0ps.value_1d(h, μ) + El0ps.conjugate_1d(h, τ) >= μ * τ - 1e-8
                else
                    @test El0ps.conjugate_1d(h, τ) < λ
                end
                @test El0ps.value(h, x) >= 0.0
                @test El0ps.value(h, z) == 0.0
                @test El0ps.value(h, x) ≈ El0ps.value(h, -x)
                @test El0ps.conjugate(h, x) >= 0.0
                @test El0ps.conjugate(h, z) == 0.0
                @test El0ps.conjugate(h, x) ≈ El0ps.conjugate(h, -x)
                @test El0ps.value(h, x) + El0ps.conjugate(h, r) >= x' * r
                @test length(El0ps.prox(h, x, η)) == length(x)
                @test El0ps.conjugate(h, x) <= Inf
                @test isa(El0ps.params_to_dict(h), OrderedDict)
            end
        end
    end
end
