function maximal_activity(row::Row{T}, j::Int, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    # compute the sup of row i (execpt the integer variable j)

    # only compute the maximal activity of any variables except j
    sup = zero(T)

    for (k, a_k) in zip(row.nzind, row.nzval)
        if k == j
            continue
        elseif a_k > zero(T)
            sup += a_k*ucol[k]
        else
            sup += a_k*lcol[k]
        end
    end
    return T(sup)
end

function minimal_activity(row::Row{T}, j::Int, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    # compute the inf of row i (execpt the integer variable j)

    # only compute the minimal activity of any variables except j
    sup = zero(T)

    for (k, a_k) in zip(row.nzind, row.nzval)
        if k == j
            continue
        elseif a_k > zero(T)
            sup += a_k*lcol[k]
        else
            sup += a_k*ucol[k]
        end
    end
    return T(sup)
end

function upperbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::Int)::T where {T}
    # perform coef strengthening for one constraints of the from a'x <= u
    row = ps.pb0.arows[i]
    a = row.nzval
    u = ps.urow[i]
    if a[j_index] > 0
        d = u - maximal_activity(row, j, ps.ucol, ps.lcol) - a[j_index]*(ps.ucol[j]-1)
        if a[j_index] >= d && d > 0
            a[j_index] = a[j_index] - d
            ps.urow[i] = u - d*ps.ucol[j]
        end
    elseif a[j_index] < 0
        d = u - maximal_activity(row, j, ps.ucol, ps.lcol) - a[j_index]*(ps.lcol[j]+1)
        if -a[j_index] >= d && d > 0
            a[j_index] = a[j_index] + d
            ps.urow[i] = u + d*ps.lcol[j]
        end
    end
    return T(a[j_index])
end

function lowerbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::Int)::T where {T}
    # perform coef strengthening for one constraints of the from l < = a'x
    row = ps.pb0.arows[i]
    a = row.nzval
    l = ps.lrow[i]
    if a[j_index] > 0
        d = -l + minimal_activity(row, j, ps.ucol, ps.lcol) + a[j_index]*(ps.lcol[j]+1)
        if a[j_index] >= d && d > 0
            a[j_index] = a[j_index] - d
            ps.lrow[i] = l - d*ps.lcol[j]
        end
    elseif a[j_index] < 0
        d = -l - minimal_activity(row, j, ps.ucol, ps.lcol) + a[j_index]*(ps.ucol[j]-1)
        if -a[j_index] >= d && d > 0
            a[j_index] = a[j_index] + d
            ps.urow[i] = u + d*ps.ucol[j]
        end
    end
    return T(a[j_index])
end

function zero_coefficient_strengthening!(ps::PresolveData{T}) where {T}
    i_index = ones(Int, ps.pb0.nvar) # keep track of index for each var fo update ps.acols

    for i in 1:ps.pb0.ncon

        row = ps.pb0.arows[i]
        lrow = ps.lrow[i]
        urow = ps.urow[i]
        if lrow > -Inf && urow < Inf #skip ranged constraints
            continue
        end

        j_index = 0 # keep track of position of variable j in row.nzind & row.nzval
        for j in row.nzind
            j_index += 1
            if ps.var_types[j] == CONTINUOUS || !ps.colflag[j]
                continue
            else
                if urow < Inf
                    new_coef = upperbound_strengthening!(ps, i, j_index, j)
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    i_index[j] += 1
                elseif lrow > -Inf
                    new_coef = lowerbound_strengthening!(ps, i, j_index, j)
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    i_index[j] += 1
                end
            end

        end
    end

    return nothing
end
