function l0screening!(
    solver::BnbSolver,
    node::BnbNode,
    F::AbstractDatafit,
    A::Matrix,
    y::Vector,
    λ::Float64,
    x::Vector, 
    w::Vector, 
    u::Vector, 
    q::Vector, 
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

    for i in findall(Sbi .| Sbb)
        if dv + λ * max(-q[i], 0.0) > ub + tolprune
            # Move i from Sbi or Sbb to S0
            Sbi[i] = false
            Sbb[i] = false
            Sb[i] = false
            S0[i] = true
            if x[i] != 0.0
                axpy!(-x[i], A[:, i], w)
                copy!(u, -gradient(F, y, w))
                x[i] = 0.0
            end
            solver.supp_pruned += 2. ^ (-depth(node))
        elseif dv + λ * max(q[i], 0.0) > ub + tolprune
            # Move i from Sbi or Sbb to S1i
            Sbi[i] = false
            Sbb[i] = false
            S1i[i] = true
            Sb[i] = false
            S1[i] = true
            solver.supp_pruned += 2. ^ (-depth(node))
        end
    end
    return nothing
end

function l1screening!(
    F::AbstractDatafit,
    A::Matrix,
    y::Vector,
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
    for i in findall(Sbi)
        vi = v[i]
        if abs(vi) + radius < τ * λ           
            # Move i from Sbi to Sb0
            if x[i] != 0.
                axpy!(-x[i], A[:, i], w)  
                copy!(u, -gradient(F, y, w))
                x[i] = 0.
            end
            Sbi[i] = false
            Sb0[i] = true
        elseif abs(vi) - radius > τ * λ        
            # Move i from Sbi to Sbb
            Sbi[i] = false
            Sbb[i] = true
        end
    end

    return nothing
end
