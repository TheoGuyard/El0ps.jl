@testset "Abstract" begin

    struct NewH <: El0ps.AbstractPerturbation end

    h = NewH()
    x = randn()
    η = rand()

    @test_throws ErrorException El0ps.value_1d(h, x)
    @test_throws ErrorException El0ps.conjugate_1d(h, x)
    @test_throws ErrorException El0ps.prox_1d(h, x, η)
end

@testset "Perturbations" begin
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
            h = test_type(test_params...)
            n = 100
            m = x = randn(n)
            z = zeros(n)
            r = randn(n)
            η = randn()
            @test isa(println(h), Nothing)
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
