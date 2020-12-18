function ceil(val::T, ϵ_int::T)::T where{T}
    if (val - floor(val)) < ϵ_int
        return T(floor(val))
    end
    return T(ceil(val))
end

function floor(val::T, ϵ_int::T)::T where{T}
    if (ceil(val) - val) < ϵ_int
        return T(ceil(val))
    end
    return T(floor(val))
end

function ensure_integer_bounds!(ps::PresolveData{T}, j::Int, ϵ_int::T=eps(T)) where{T}
    # Column was already removed or the variable is continous.
    ps.colflag[j] || (ps.var_types[j] != CONTINUOUS) || return nothing

    if (ps.var_types[j] == BINARY)
        ps.lcol[j] = T(max(0, ceil(ps.lcol[j], ϵ_int)))
         ps.ucol[j] = T(min(1, floor(ps.ucol[j], ϵ_int)))
    elseif (ps.var_types[j] == GENERAL_INTEGER)
        ps.lcol[j], ps.ucol[j] = ceil(ps.lcol[j], ϵ_int), floor(ps.ucol[j], ϵ_int)
        if (ps.lcol[j] >= T(0)) && (ps.ucol[j] <= T(1))
            ps.var_types[j] == BINARY
        end
    end

    # Only check the infeasibility for integer variables.
    if (ps.lcol[j] > ps.ucol[j])
        ps.status = PRIMAL_INFEASIBLE
    end
    return nothing
end
