Base.@kwdef mutable struct PresolveOptions{T}
    Level::Int = 1  # Presolve level

    # Numerical tolerances
    PTol::T = eps(T)  # Primal tolerance
    DTol::T = eps(T)  # Dual tolerance
    ITol::T = T(1 // 10_000)  # Integrality tolerance (1e-5)
end
