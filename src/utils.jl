function bisection(f::Function, a::Float64, b::Float64; ϵ::Float64=1e-8, maxit::Int=100)
    it = 1
    while it <= maxit
        c = (a + b) / 2.0
        fa = f(a)
        fc = f(c)
        (fc ≈ 0.0) && return c
        (b - a < ϵ) && return c
        it += 1
        if fc * fa >= 0.
            a = c
        else
            b = c
        end
    end
    return c
end

function approximate_τ(h::AbstractPenalty)
    a = 0.0
    b = 1.0
    while conjugate_1d(h, b) < 1.0
        b *= 2.0
    end
    return bisection(v -> conjugate_1d(h, v) - 1.0, a, b)
end
