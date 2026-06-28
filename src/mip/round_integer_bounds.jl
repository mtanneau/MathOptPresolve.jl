@doc raw"""
    RoundSingleIntegerBound <: AbstractRule

Round the upper and lower bounds of single integer variable to be integers,
i.e. lcol[j] = ceil(lcol[j]) and ucol[j] = floor(lcol[j]) for general integers
and lcol[j] = max(0, ceil(lcol[j])) and ucol[j] = min(1, floor(lcol[j])) for
binary variables.
"""

struct RoundSingleIntegerBound <: AbstractRule
    j::Int
end

function apply!(
    ps::PresolveData{T},
    r::RoundSingleIntegerBound,
    config::PresolveOptions{T}
) where {T}
    j = r.j
    ϵ_int = config.IntegerTolerance

    # Column was already removed or the variable is continous.
    ps.colflag[j] || return nothing
    (ps.var_types[j] == CONTINUOUS) && return nothing

    if (ps.var_types[j] == BINARY) # Variable type is binary.
        ps.lcol[j] = T(max(0, approx_ceil(ps.lcol[j], ϵ_int)))
        ps.ucol[j] = T(min(1, approx_floor(ps.ucol[j], ϵ_int)))
    elseif (ps.var_types[j] == GENERAL_INTEGER) # Variable type is integer.
        ps.lcol[j] = approx_ceil(ps.lcol[j], ϵ_int)
        ps.ucol[j] = approx_floor(ps.ucol[j], ϵ_int)

        if (ps.lcol[j] >= T(0)) && (ps.ucol[j] <= T(1)) # It is binary.
            ps.var_types[j] = BINARY
        end
    end

    return nothing
end


"""
    RoundIntegerBounds <: AbstractRule

Round the upper and lower bounds of all integer variables to be integers.
"""

struct RoundIntegerBounds <: AbstractRule end

function apply!(
    ps::PresolveData{T},
    ::RoundIntegerBounds,
    config::PresolveOptions{T}
) where {T}
    # The problem is LP.
    ps.pb0.is_continuous && return nothing

    for j in 1:ps.pb0.nvar
        apply!(ps, RoundSingleIntegerBound(j), config)
    end
    return nothing
end
