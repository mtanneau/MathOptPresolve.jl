function maximal_activity(row::RowOrCol{T}, j::Int, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    # compute the sup of row i (execpt the integer variable j)

    # only compute the maximal activity of any variables except j
    nonzero_index = [row.nzind[i] for i in 1:length(row.nzind) if row.nzind[i] != j]
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

function single_row_strengthening(row::RowOrCol{T}, b::T, j::Int, ucol::Vector{T}, lcol::Vector{T}) where {T}
    # perform coef strengthening for one constraints of the from a'x < b
    # and return new a' and b'
    flag = false # wether there is an update
    i = findfirst(isequal(j), row.nzind) # index for integer variable j
    a = deepcopy(row.nzval)
    if a[i] > 0
        d = b - maximal_activity(row, j, ucol, lcol) - a[i]*(ucol[j]-1)
        if a[i] >= d && d > 0
            a[i] = a[i] - d
            b = b - d*ucol[j]
            flag = true
        end
    elseif a[i] < 0
        d = b - maximal_activity(row, j, ucol, lcol) - a[i]*(lcol[j]+1)
        if -a[i] >= d && d > 0
            a[i] = a[i] + d
            b = b + d*lcol[j]
            flag = true
        end
    end
    return a, b, a[i], flag
end

function coefficient_strengthening!(ps::PresolveData{T}, j::Int, inf::T = T(1e30)) where {T}
    # perform coefficient strengthening on integer variable j
    ps.colflag[j] || (ps.var_types[j] != CONTINUOUS) || return nothing

    for i in 1:length(ps.pb0.arows)
        row = ps.pb0.arows[i]
        lrow = ps.lrow[i]
        urow = ps.urow[i]

        if lrow > -inf && urow < inf
            continue #skipping ranged constraints
        elseif urow < inf
            a, b, new_coef, updated = single_row_strengthening(row, urow, j, ps.ucol, ps.lcol)
            if updated
                row.nzval = a
                ps.urow[i] = b

                # update collumns of A in problem data
                k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                ps.pb0.acols[j].nzval[k] = new_coef
            end

        elseif lrow > -inf
            r = deepcopy(row)
            r.nzval = -r.nzval
            a, b, new_coef, updated = single_row_strengthening(r, -lrow, j, ps.ucol, ps.lcol)
            if updated
                row.nzval = -a
                ps.lrow[i] = -b

                # update collumns of A in problem data
                k = findfirst(isequal(i), ps.pb0.acols[j].nzind)
                ps.pb0.acols[j].nzval[k] = -new_coef
            end
        end
    end

    return nothing
end
