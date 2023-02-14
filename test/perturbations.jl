@testset "Perturbations" begin
    candidates = [
        (Zero, ()),
        (Bigm, (1.)),
        (L1norm, (1.)),
        (L2norm, (1.)),
        (L1L2norm, (1., 1.)),
        (BigmL1norm, (1., 1.)),
        (BigmL2norm, (1., 1.)),
    ]
    for (test_type, test_params) in candidates
        @testset "$test_type utilities" begin
            h = test_type(test_params...)
            n = 100
            m = 
            x = randn(n)
            z = zeros(n)
            r = randn(n)
            η = randn()
            @test El0ps.value(h, x) >= 0.
            @test El0ps.value(h, z) == 0.
            @test El0ps.value(h, x) ≈ El0ps.value(h, -x)
            @test El0ps.conjugate(h, x) >= 0.
            @test El0ps.conjugate(h, z) == 0.
            @test El0ps.conjugate(h, x) ≈ El0ps.conjugate(h, -x)
            @test El0ps.value(h, x) + El0ps.conjugate(h, r) >= x' * r
            @test length(El0ps.prox(h, x, η)) == length(x)
            @test El0ps.conjugate(h, x) <= Inf
            @test isa(El0ps.params_to_dict(h), OrderedDict)
        end
    end
end
