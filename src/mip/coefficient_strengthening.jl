function maximal_activity(row::Row{T}, j::Int, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    # compute the sup of row i (execpt the integer variable j)

    # only compute the maximal activity of any variables except j
    sup = zero(T)

    for (k, a_k) in zip(row.nzind, row.nzval)
        if k == j
            continue
        elseif a_k > zero(T)
            sup += a_k*ucol[j]
        else
            sup += a_k*lcol[j]
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
            sup += a_k*lcol[j]
        else
            sup += a_k*ucol[j]
        end
    end
    return T(sup)
end

function upperbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::In) where {T}
    # perform coef strengthening for one constraints of the from a'x <= u
    row = ps.pb0.arows[i]
    a = row.nzval
    u = ps.urow[i]
    if a[j_index] > 0
        d = u - maximal_activity(row, j, ps.ucol, ps.lcol) - a[j_index]*(ps.ucol[j]-1)
        if a[j_index] >= d && d > 0
            a[j_index] = a[j_index] - d
            ps.urow[i] = u - d*ucol[j]
        end
    elseif a[j_index] < 0
        d = u - maximal_activity(row, j, ps.ucol, ps.lcol) - a[j_index]*(ps.lcol[j]+1)
        if -a[j_index] >= d && d > 0
            a[j_index] = a[j_index] + d
            ps.urow[i] = u + d*lcol[j]
        end
    end
    return nothing
end

function lowerbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::In) where {T}
    # perform coef strengthening for one constraints of the from l < = a'x
    row = ps.pb0.arows[i]
    a = row.nzval
    l = ps.urow[i]
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
    return nothing
end

function coefficient_strengthening!(ps::PresolveData{T}, i::Int) where {T}
    # perform coefficient strengthening on row i
    row = ps.pb0.arows[i]
    lrow = ps.lrow[i]
    urow = ps.urow[i]

    if lrow > -Inf && urow < Inf #skip ranged constraints
        return nothing
    end

    for j in ps.pb0.nvar
        if ps.var_types[j] == CONTINUOUS || !ps.colflag[j]
            continue
        else
            j_index = findfirst(isequal(j), row.nzind) # index of var j in sparse row
            if j_index == nothing # zero coefficent variable
                continue
            elseif urow < Inf
                upperbound_strengthening!(ps, i, j_index, j)
            elseif lrow > -Inf
                lowerbound_strengthening!(ps, i, j_index, j)
            end
        end
    end

    return nothing
end
