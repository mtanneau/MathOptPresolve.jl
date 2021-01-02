function bound_strengthening!(ps::PresolveData{T}, j::Int, ϵ::T=eps(T), ϵ_int::T=eps(T)) where{T}
    # Column was already removed or the variable is continous.
    ps.colflag[j] || return nothing
    (ps.var_types[j] == CONTINUOUS) && return nothing

    for i in 1:ps.pb0.ncon
        row = ps.pb0.arows[i]
        j_ind = findfirst(isequal(j), row.nzind)
        if (j_ind == nothing)
            row_j = 0
        else
            row_j = row.nzval[j_ind]
        end
        if (abs(row_j) <= ϵ) # a_{i,j} is 0.
            continue
        end
        lrow = ps.lrow[i]
        urow = ps.urow[i]
        upper = calc_upper_bound_except_one(row, ps.lcol, ps.ucol, j)
        lower = calc_upper_bound_except_one(row, ps.lcol, ps.ucol, j)
        if (row_j > T(0))
            if (urow - lower != NaN)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((urow - lower) / row_j, ϵ_int)))
            end
            if (lrow - upper != NaN)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((lrow - upper) / row_j, ϵ_int)))
            end
        else
            if (lrow - upper != NaN)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((urow - lower) / row_j, ϵ_int)))
            end
            if (urow - lower != NaN)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((lrow - upper) / row_j, ϵ_int)))
            end
        end
    end
    return nothing
end

function calc_upper_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        (val > 0) ? bound += val * ucol[ind] : bound += val * lcol[ind]
    end
    return bound
end

function calc_lower_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        (val > 0) ? bound += val * lcol[ind] : bound += val * ucol[ind]
    end
    return bound
end
