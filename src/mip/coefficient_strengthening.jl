function maximal_activity(row::Row{T}, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    sup = zero(T)

    for (j, aij) in zip(row.nzind, row.nzval)
        if aij > zero(T)
            sup += aij*ucol[j]
        else
            sup += aij*lcol[j]
        end
    end
    return T(sup)
end

function minimal_activity(row::Row{T}, ucol::Vector{T}, lcol::Vector{T})::T where {T}
    inf = zero(T)

    for (j, aij) in zip(row.nzind, row.nzval)
        if aij > zero(T)
            inf += aij*lcol[j]
        else
            inf += aij*ucol[j]
        end
    end
    return T(inf)
end

function upperbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::Int, max_act) where {T}
    # perform coef strengthening for one constraints of the from a'x <= u
    row = ps.pb0.arows[i]
    a = row.nzval
    u = ps.urow[i]
    if a[j_index] > 0
        d = u - max_act - a[j_index]*(ps.ucol[j]-1)
        if a[j_index] >= d && d > 0
            a[j_index] = a[j_index] - d
            ps.urow[i] = u - d*ps.ucol[j]
        end
    elseif a[j_index] < 0
        d = u - max_act - a[j_index]*(ps.lcol[j]+1)
        if -a[j_index] >= d && d > 0
            a[j_index] = a[j_index] + d
            ps.urow[i] = u + d*ps.lcol[j]
        end
    end
    return T(a[j_index])
end

function lowerbound_strengthening!(ps::PresolveData{T}, i::Int, j_index::Int, j::Int, min_act) where {T}
    # perform coef strengthening for one constraints of the from l < = a'x
    row = ps.pb0.arows[i]
    a = row.nzval
    l = ps.lrow[i]
    if a[j_index] > 0
        d = -l + min_act + a[j_index]*(ps.lcol[j]+1)
        if a[j_index] >= d && d > 0
            a[j_index] = a[j_index] - d
            ps.lrow[i] = l - d*ps.lcol[j]
        end
    elseif a[j_index] < 0
        d = -l - min_act + a[j_index]*(ps.ucol[j]-1)
        if -a[j_index] >= d && d > 0
            a[j_index] = a[j_index] + d
            ps.urow[i] = u + d*ps.ucol[j]
        end
    end
    return T(a[j_index])
end

function zero_coefficient_strengthening!(ps::PresolveData{T}) where {T}
    # perform coefficient stregthening but if there is a coefficient is reduced to 0
    # it is still kept in the sparse matrix (will be removed later)

    # keep track of index for each var fo update ps.acols
    # use this to find which index of ps.pb0.acols[i].nzval to update

    i_index = zeros(Int, ps.pb0.nvar)
    for i in 1:ps.pb0.ncon

        lrow = ps.lrow[i]
        urow = ps.urow[i]
        if lrow > -Inf && urow < Inf #skip ranged constraints
            continue
        end

        row = ps.pb0.arows[i]
        sup = maximal_activity(row, ps.ucol, ps.lcol)
        inf = minimal_activity(row, ps.ucol, ps.lcol)

        j_index = 0 # keep track of index of variable j in row.nzind & row.nzval
        for j in row.nzind
            j_index += 1
            i_index[j] += 1
            if ps.var_types[j] == CONTINUOUS || !ps.colflag[j]
                continue
            else
                coef = row.nzval[j_index]
                if urow < Inf
                    if coef > 0
                        max_act = sup - coef * ps.ucol[j] # maximal activity of every variables except j
                    else
                        max_act = sup - coef * ps.lcol[j]
                    end
                    new_coef = upperbound_strengthening!(ps, i, j_index, j, max_act)
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    #update sup
                    if coef > 0
                        sup = sup - (coef - new_coef) * ps.ucol[j]
                    else
                        sup = sup - (coef - new_coef) * ps.lcol[j]
                    end
                elseif lrow > -Inf
                    if coef > 0
                        min_act = inf - coef * ps.lcol[j] # minimal activity of every variables except j
                    else
                        min_act = inf - coef * ps.ucol[j]
                    end
                    new_coef = lowerbound_strengthening!(ps, i, j_index, j, min_act)
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    #update inf
                    if coef > 0
                        inf = inf - (coef - new_coef) * ps.lcol[j]
                    else
                        inf = inf - (coef - new_coef) * ps.ucol[j]
                    end
                end
            end

        end
    end
    return nothing
end
