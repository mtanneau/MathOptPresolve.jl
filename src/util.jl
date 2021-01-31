# Positive and negative part of a number
pos_part(x::T) where {T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where {T} = x >= zero(T) ? zero(T) : -x

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
