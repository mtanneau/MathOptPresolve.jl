# Positive and negative part of a number
pos_part(x::T) where {T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where {T} = x >= zero(T) ? zero(T) : -x
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

function maximal_activity(ps::PresolveData{T}, i::Int)::T where {T}
    sup = zero(T)

    row = ps.pb0.arows[i]
    for (j, aij) in zip(row.nzind, row.nzval)
        ps.colflag[j] || continue
        if aij > zero(T)
            sup += aij * ps.ucol[j]
        elseif aij < zero(T)
            sup += aij * ps.lcol[j]
        end
        isfinite(sup) || return sup
    end
    return sup
end

function minimal_activity(ps::PresolveData{T}, i::Int)::T where {T}
    inf = zero(T)

    row = ps.pb0.arows[i]
    for (j, aij) in zip(row.nzind, row.nzval)
        ps.colflag[j] || continue
        if aij > zero(T)
            inf += aij * ps.lcol[j]
        elseif aij < zero(T)
            inf += aij * ps.ucol[j]
        end
        isfinite(inf) || return inf
    end
    return inf
end
