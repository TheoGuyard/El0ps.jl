@testset "Utils" begin

    @testset "Bisection" begin
        f = x -> 3.0 * x^3 - 4.0 * x^2 + 2.0 * x - 1.0
        ϵ = 1e-8
        c = El0ps.bisection(f, 0.0, 2.0, ϵ = ϵ)
        @assert abs(c - 1.0) < ϵ
        @assert f(c) ≈ 0.0
    end

    @testset "Approximate τ" begin
        candidates = [
            (El0ps.Bigm, (1.1)),
            (El0ps.L1norm, (1.1)),
            (El0ps.L2norm, (1.1)),
            (El0ps.L1L2norm, (1.1, 1.2)),
            (El0ps.BigmL1norm, (1.1, 1.2)),
            (El0ps.BigmL2norm, (1.1, 1.2)),
            (El0ps.NegLogSymtri, (1.1, 1.2)),
        ]
        for (test_type, test_params) in candidates
            @testset "With $test_type" begin
                h = test_type(test_params...)
                ϵ = 1e-8
                τ1 = El0ps.compute_τ(h)
                τ2 = El0ps.approximate_τ(h, ϵ = ϵ)
                @assert abs(τ1 - τ2) < ϵ
            end
        end
    end


end

function bisection(f::Function, a::Float64, b::Float64; ϵ::Float64 = 1e-8, maxit::Int = 100)

end

function approximate_τ(h::AbstractPenalty)
    a = 0.0
    b = 1.0
    while conjugate_1d(h, b) < 1.0
        b *= 2.0
    end
    return bisection(v -> conjugate_1d(h, v) - 1.0, a, b)
end
