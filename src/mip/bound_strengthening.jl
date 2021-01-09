@doc raw"""

The behavior of bound_strengthening!() depends on the ordering of the variables.

Let the constraints in the j-th row be
lrow ⩽ row_S x_S + row_j x_j ⩽ urow

Denote the U(row_S x_S) be the upper bound of row_S x_S and
L(row_S x_S) be the lower bound of row_S x_S.

Then, lrow - U(row_S x_S) ⩽ row_j x_j ⩽ urow - L(row_S x_S).

Thus, if row_j > 0,
x_j ⩽ floor((urow - L(row_S x_S)) / row_j) and
x_j ⩾ ceil((lrow - U(row_S x_S) x_j) / row_j).

If row_j < 0,
x_j ⩾ ceil((urow - L(row_S x_S)) / row_j) and
x_j ⩽ floor((lrow - U(row_S x_S) x_j) / row_j).
"""

function bound_strengthening!(ps::PresolveData{T}, j::Int, ϵ_int::T=eps(T)) where{T}
    # Column was already removed or the variable is continous.
    ps.colflag[j] || return nothing
    (ps.var_types[j] == CONTINUOUS) && return nothing
    col = ps.pb0.acols[j]
    for (i, row_j) in zip(col.nzind, col.nzval)
        # Row was already removed.
        ps.rowflag[i] || continue
        row = ps.pb0.arows[i]

        lrow = ps.lrow[i]
        urow = ps.urow[i]
        upper = calc_upper_bound_except_one(row, ps.lcol, ps.ucol, j)
        lower = calc_lower_bound_except_one(row, ps.lcol, ps.ucol, j)
        if (row_j > T(0))
            if isfinite(urow) && isfinite(lower)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((urow - lower) / row_j, ϵ_int)))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((lrow - upper) / row_j, ϵ_int)))
            end
        else
            if isfinite(urow) && isfinite(lower)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((urow - lower) / row_j, ϵ_int)))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((lrow - upper) / row_j, ϵ_int)))
            end
        end
    end
    return nothing
end
