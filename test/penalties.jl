@testset "Penalties" begin
    candidates = [
        (Bigm, (1.)),
        (L1norm, (1.)),
        (L2norm, (1.)),
        (L1L2norm, (1., 1.)),
        (BigmL1norm, (1., 1.)),
        (BigmL2norm, (1., 1.)),
    ]

    for (test_type, test_params) in candidates
        G = test_type(test_params...)
        @testset "$G" begin
            n = 100
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
        end
    end
end