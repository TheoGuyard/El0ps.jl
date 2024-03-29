function l0screening!(
    f::AbstractDatafit,
    A::Matrix,
    λ::Float64,
    x::Vector,
    w::Vector,
    u::Vector,
    p::Vector,
    ub::Float64,
    dv::Float64,
    tolprune::Float64,
    S0::BitArray,
    S1::BitArray,
    Sb::BitArray,
    Sbi::BitArray,
    Sbb::BitArray,
    S1i::BitArray,
)

    idx_to_test = findall(@. Sb & (!isnan(p)))
    for i in idx_to_test
        if dv + λ * max(-p[i], 0.0) > ub + tolprune
            # Move i from Sbi or Sbb to S0
            Sbi[i] = false
            Sbb[i] = false
            Sb[i] = false
            S0[i] = true
            if x[i] != 0.0
                axpy!(-x[i], A[:, i], w)
                copy!(u, -gradient(f, w))
                x[i] = 0.0
            end
        elseif dv + λ * max(p[i], 0.0) > ub + tolprune
            # Move i from Sbi or Sbb to S1i
            Sbi[i] = false
            Sbb[i] = false
            S1i[i] = true
            Sb[i] = false
            S1[i] = true
        end
    end
    return nothing
end

function l1screening!(
    f::AbstractDatafit,
    A::Matrix,
    λ::Float64,
    x::Vector,
    w::Vector,
    u::Vector,
    v::Vector,
    α::Float64,
    τ::Float64,
    gap::Float64,
    Sb0::BitArray,
    Sbi::BitArray,
    Sbb::BitArray,
)

    radius = sqrt(2.0 * gap * α)
    idx_to_test = findall(@. Sbi & !isnan(v))
    for i in idx_to_test
        vi = v[i]
        if abs(vi) + radius < λ * τ
            # Move i from Sbi to Sb0
            if x[i] != 0.0
                axpy!(-x[i], A[:, i], w)
                copy!(u, -gradient(f, w))
                x[i] = 0.0
            end
            Sbi[i] = false
            Sb0[i] = true
        elseif abs(vi) - radius > λ * τ
            # Move i from Sbi to Sbb
            Sbi[i] = false
            Sbb[i] = true
        end
    end

    return nothing
end
