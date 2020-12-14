Base.@kwdef mutable struct PresolveOptions{T}
    Level::Int = 1  # Presolve level

    # Numerical tolerances
    PrimalTolerance::T = eps(T)  # Primal tolerance
    DualTolerance::T = eps(T)  # Dual tolerance
    IntegerTolerance::T = T(1 // 10_000)  # Integrality tolerance (1e-5)
end
