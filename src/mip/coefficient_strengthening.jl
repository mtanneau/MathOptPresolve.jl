function maximal_activity(row::RowOrCol{T}, j::Int, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    # compute the sup of row i (execpt the integer variable j)

    # only compute the maximal activity of any variables except j
    nonzero_index = [i for i in 1:length(row.nzind) if row.nzind[i] != j]
    nonzero_value = row.nzval[nonzero_index]

    pos_coef = deepcopy(nonzero_value)
    neg_coef = deepcopy(nonzero_value)

    for i in 1:length(nonzero_value)
        if nonzero_value[i] < 0
            pos_coef[i] = 0
        else
            neg_coef[i] = 0
        end
    end
    sup = ucol[nonzero_index]'pos_coef + lcol[nonzero_index]'neg_coef
    return T(sup)
end

function single_row_strengthening(row::RowOrCol{T}, b::T, m::Int, j::Int, ucol::Vector{T}, lcol::Vector{T}) where {T}
    # perform coef strengthening for one constraints of the from a'x < b
    # and return new a' and b'
    flag = false # wether there is an update
    a = deepcopy(row.nzval)
    if a[m] > 0
        d = b - maximal_activity(row, j, ucol, lcol) - a[m]*(ucol[j]-1)
        if a[m] >= d && d > 0
            a[m] = a[m] - d
            b = b - d*ucol[j]
            flag = true
        end
    elseif a[m] < 0
        d = b - maximal_activity(row, j, ucol, lcol) - a[m]*(lcol[j]+1)
        if -a[m] >= d && d > 0
            a[m] = a[m] + d
            b = b + d*lcol[j]
            flag = true
        end
    end
    return a, b, a[m], flag
end

function coefficient_strengthening!(ps::PresolveData{T}, j::Int) where {T}
    # perform coefficient strengthening on integer variable j
    (ps.var_types[j] != CONTINUOUS) || return nothing

    for i in 1:length(ps.pb0.arows)
        row = ps.pb0.arows[i]
        lrow = ps.lrow[i]
        urow = ps.urow[i]

        m = findfirst(isequal(j), row.nzind)
        if m == nothing || (lrow > -Inf && urow < Inf)
            continue #skipping ranged constraints and
        elseif urow < Inf
            a, b, new_coef, updated = single_row_strengthening(row, urow, m, j, ps.ucol, ps.lcol)
            if updated
                row.nzval = a
                ps.urow[i] = b

                # update collumns of A in problem data
                # and update the sparse matrix in the case where the new coefficient is close to 0
                if abs(new_coef) <= eps(T)
                    deleteat!(row.nzind, m)
                    deleteat!(row.nzval, m)
                    k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                    deleteat!(ps.pb0.acols[j].nzind, k)
                    deleteat!(ps.pb0.acols[j].nzval, k)
                else
                    k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                    ps.pb0.acols[j].nzval[k] = new_coef
                end
            end
        elseif lrow > -Inf
            r = deepcopy(row)
            r.nzval = -r.nzval
            a, b, new_coef, updated = single_row_strengthening(r, -lrow, m, j, ps.ucol, ps.lcol)
            if updated
                row.nzval = -a
                ps.lrow[i] = -b

                # update collumns of A in problem data
                if abs(new_coef) <= eps(T)
                    deleteat!(row.nzind, m)
                    deleteat!(row.nzval, m)
                    k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                    deleteat!(ps.pb0.acols[j].nzind, k)
                    deleteat!(ps.pb0.acols[j].nzval, k)
                else
                    k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                    ps.pb0.acols[j].nzval[k] = -new_coef
                end
            end
        end
    end
    return nothing
end
