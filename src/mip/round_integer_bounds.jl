function round_integer_bounds!(ps::PresolveData{T}, j::Int, ϵ_int::T=eps(T)) where{T}
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
