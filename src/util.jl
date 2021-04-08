# Positive and negative part of a number
pos_part(x::T) where {T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where {T} = x >= zero(T) ? zero(T) : -x

# Maximal activity of a constraint
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

# Minimal activity of a constraint
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

"""
In some presolve procedures, when we process rows, we need to update
columns accordingly. This function performs synchronization between
rows and columns after a presolve method is applied on every rows.
"""

function sync_columns_to_rows!(ps::PresolveData{T}) where {T}
    # keep track of how far to update columns
    i_index = zeros(Int, ps.pb0.nvar)
    for i in 1:ps.pb0.ncon
        nzind_i = ps.pb0.arows[i].nzind
        for j in nzind_i
            i_index[j] += 1
        end

        for (k, j) in enumerate(nzind_i)
            ps.pb0.acols[j].nzval[i_index[j]] = ps.pb0.arows[i].nzval[k]
        end
    end
    return nothing
end
