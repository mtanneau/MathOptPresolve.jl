"""
Perform Chavatal-Gomory strengthening of inequalities
"""
function cg_strengthening!(ps::PresolveData{T}) where {T}
    ps.pb0.is_continuous && return nothing

    i_index = zeros(Int, ps.pb0.nvar)
    for i in 1:ps.pb0.ncon
        nonzero_index = ps.pb0.arows[i].nzind
        i_index[nonzero_index] .+= 1
        ps.rowflag[i] || continue

        if isfinite(ps.urow[i]) && isfinite(ps.lrow[i]) ||
            any(ps.var_types[nonzero_index] .== CONTINUOUS)
            # skip ranged constraints
            # and any constraint that involves continuous variables
            continue
        elseif isfinite(ps.lrow[i])
            lowerbound_cg!(ps, i, i_index)
        elseif isfinite(ps.urow[i])
            upperbound_cg!(ps, i, i_index)
        end
    end
    return nothing
end

function lowerbound_cg!(ps::PresolveData{T}, i::Int, i_index::Vector{Int}, search_range::Int = 5) where {T}
    # perform CG strengthening inequality for a single constraint of the
    # from a'x ⩾ b
    l = ps.lrow[i]
    row = ps.pb0.arows[i].nzval
    nonzero_index = ps.pb0.arows[i].nzind
    row_pos_index = []
    row_neg_index = []
    col_pos_index = []
    col_neg_index = []
    for j in enumerate(nonzero_index)
        if row[j[1]] > 0 && isfinite(ps.lcol[j[2]])
            l -= row[j[1]] * ps.lcol[j[2]]
            append!(row_pos_index, j[1])
            append!(col_pos_index, j[2])
        elseif row[j[1]] < 0 && isfinite(ps.ucol[j[2]])
            l -= row[j[1]] * ps.ucol[j[2]]
            append!(row_neg_index, j[1])
            append!(col_neg_index, j[2])
        else
            return nothing
        end
    end

    # building set for heuristicly finding s
    a = abs.(row)
    a_max = maximum(a)
    a_min = minimum(a)

    set = Set(T(1))
    for t in 1:search_range
        union!(set, [t/a_max, t/a_min, (2*t-1)/(2*a_min)])
    end

    # perform search for s and update
    for s in set
        if any(a*ceil(l*s) < ceil.(s*a)*l)
            continue
        elseif any(a*ceil(l*s) > ceil.(s*a)*l)
            #update row
            row[row_pos_index] = ceil.(s*a[row_pos_index])
            row[row_neg_index] = -ceil.(s*a[row_neg_index])

            #update column
            for j in enumerate(nonzero_index)
                ps.pb0.acols[j[2]].nzval[i_index[j[2]]] = row[j[1]]
            end

            #update bound
            ps.lrow[i] = ceil(s*l) +
            ceil.(s*a[row_pos_index])'ps.lcol[col_pos_index] -
            ceil.(s*a[row_neg_index])'ps.ucol[col_neg_index]
            break
        end
    end
    return nothing
end

function upperbound_cg!(ps::PresolveData{T}, i::Int, i_index::Vector{Int}, search_range::Int = 5) where {T}
    # perform CG strengthening inequality for a single constraint of the
    # from a'x ⩽ b
    u = ps.urow[i]
    row = ps.pb0.arows[i].nzval
    nonzero_index = ps.pb0.arows[i].nzind
    row_pos_index = []
    row_neg_index = []
    col_pos_index = []
    col_neg_index = []
    for j in enumerate(nonzero_index)
        if row[j[1]] > 0 && isfinite(ps.ucol[j[2]])
            u -= row[j[1]] * ps.ucol[j[2]]
            append!(row_pos_index, j[1])
            append!(col_pos_index, j[2])
        elseif row[j[1]] < 0 && isfinite(ps.lcol[j[2]])
            u -= row[j[1]] * ps.lcol[j[2]]
            append!(row_neg_index, j[1])
            append!(col_neg_index, j[2])
        else
            return nothing
        end
    end

    # building set for heuristicly finding s
    a = abs.(row) # assume row is nonzero
    a_max = maximum(a)
    a_min = minimum(a)

    set = Set(T(1))
    for t in 1:search_range
        union!(set, [t/a_max, t/a_min, (2*t-1)/(2*a_min)])
    end

    # perform search for s and update
    for s in set
        if any(a*floor(u*s) > floor.(s*a)*u)
            continue
        elseif any(a*floor(u*s) < floor.(s*a)*u)
            #update row
            row[row_pos_index] = ceil.(s*a[row_pos_index])
            row[row_neg_index] = floor.(-s*a[row_neg_index])

            #update column
            for j in enumerate(nonzero_index)
                ps.pb0.acols[j[2]].nzval[i_index[j[2]]] = row[j[1]]
            end

            #update bound
            ps.urow[i] = floor(s*u) +
            ceil.(s*a[row_pos_index])'ps.ucol[col_pos_index] +
            floor.(-s*a[row_neg_index])'ps.lcol[col_neg_index]
            break
        end
    end
    return nothing
end
