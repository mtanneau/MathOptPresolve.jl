# Allow some numerical errors to avoid the situation
# like ceiling 1e-10 to 1.
function approx_ceil(val::T, ϵ_int::T)::T where{T}
    if (val - floor(val)) < ϵ_int
        return T(floor(val))
    end
    return T(ceil(val))
end

# Allow some numerical errors to avoid the situation
# like flooring 1 - 1e-10 to 0.
function approx_floor(val::T, ϵ_int::T)::T where{T}
    if (ceil(val) - val) < ϵ_int
        return T(ceil(val))
    end
    return T(floor(val))
end

function round_integer_bounds!(ps::PresolveData{T}, j::Int, ϵ_int::T=eps(T)) where{T}
    # Column was already removed or the variable is continous.
    ps.colflag[j] || (ps.var_types[j] != CONTINUOUS) || return nothing

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
