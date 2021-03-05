"""
Perform Chavatal-Gomory strengthening of inequalities

We apply the Chvatal-Goromory procedure to constraints that
satisfy the following criterions:
1. All variables with nonzero coefficient must be integral
2. The constraints are not ranged
3. If the constraints is of the form ax ⩾ b, then every
   positive coefficient variables must have finite lowerbound
   negative coefficient variables must have finite upperbound

   If the constraints is of the from ax ⩽ b, then every
   positive coefficient variables must have finite upperbound
   negative coefficient variables must have finite lowerbound

For example, let aᵢx ⩾ bᵢ is a constraint where aᵢ ⩾ 0, and x
is integer variable. If there exist s > 0 such that
⌈aᵢⱼ s⌉ bᵢ / ⌈bᵢ s⌉ ⩽ aᵢⱼ
and there is k such that
⌈aᵢₖ s⌉ bᵢ / ⌈bᵢ s⌉ < aᵢₖ (*)
then we will update aᵢx ⩾ bᵢ into ⌈aᵢ s⌉ x ⩾ ⌈bᵢ s⌉

We find such a number s heuristicaly i.e, we try different
values from the following set
{1, t/a_max, t/a_min, (2t-1)/(2a_min)| t = 1,...,search_range}
where a_max = sup(|aᵢ|), a_min = inf(|aᵢ|) and search_range
is set to be 5 by default

"""
function cg_strengthening!(ps::PresolveData{T}, search_range::Int = 5) where {T}
    ps.pb0.is_continuous && return nothing

    # building a heuristic set S, for every constraints, we will try every every
    # value in this set to see if any of them satisfies condition (*)

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
            row = ps.pb0.arows[i]
            cg_inequality!(ps, row, ps.lcol, ps.ucol, ps.lrow[i], i, i_index, ceil)
        elseif isfinite(ps.urow[i])
            row = ps.pb0.arows[i]
            cg_inequality!(ps, row, ps.ucol, ps.lcol, ps.urow[i], i, i_index, floor)
        end
    end
    return nothing
end

function build_heuristic_set(a::Vector{T}, search_range::Int = 5) where {T}
    # building a heuristic set, for every constraints, we will try every every
    # value in this set to see if any of them satisfies condition (*)
    a_max = maximum(a)
    a_min = minimum(a)

    set = [T(1)]
    for t in 1:search_range
        append!(set, [t/a_max, t/a_min, (2*t-1)/(2*a_min)])
    end
    return unique(set)
end

function cg_inequality!(ps::PresolveData{T},
    row::Row{T}, var_bound_1::Vector{T}, var_bound_2::Vector{T},
    b::T, i::Int, i_index::Vector{Int}, round_function) where {T}
    # perform CG strengthening inequality for a single constraint

    # if var_bound_1 is variable lower bound, var_bound_2 is variable upper bound
    # and round_function is ceil(), we will perform CG strengthening on a constraint
    # of the form ax ⩾ b. Otherwise, we will perfom GC strengthening on a constraint
    # of the form ax ⩽ b.

    r = row.nzval
    nonzero_index = row.nzind
    row_pos_index = []
    row_neg_index = []
    col_pos_index = []
    col_neg_index = []
    for j in enumerate(nonzero_index)
        if r[j[1]] > 0 && isfinite(var_bound_1[j[2]])
            b -= r[j[1]] * var_bound_1[j[2]]
            append!(row_pos_index, j[1])
            append!(col_pos_index, j[2])
        elseif r[j[1]] < 0 && isfinite(var_bound_2[j[2]])
            b -= r[j[1]] * var_bound_2[j[2]]
            append!(row_neg_index, j[1])
            append!(col_neg_index, j[2])
        else
            return nothing
        end
    end

    # building set for heuristicly finding s
    a = abs.(r)
    set = build_heuristic_set(a)

    # perform search for s and update
    if round_function == ceil
        sign = 1
    else
        sign = -1
    end
    for s in set
        found_s = false
        for val in sign*a
            if val*round_function(b*s) < round_function(s*val)*b
                found_s = false
                break
            elseif val*round_function(b*s) > round_function(s*val)*b
                found_s = true
            end
        end

        if found_s
            # update row
            r[row_pos_index] = ceil.(s*a[row_pos_index])
            r[row_neg_index] = -ceil.(s*a[row_neg_index])

            #update column
            for j in enumerate(nonzero_index)
                ps.pb0.acols[j[2]].nzval[i_index[j[2]]] = r[j[1]]
            end

            #update bound
            new_bound = round_function(s*b) +
            ceil.(s*a[row_pos_index])'var_bound_1[col_pos_index] -
            ceil.(s*a[row_neg_index])'var_bound_2[col_neg_index]

            if sign == 1
                ps.lrow[i] = new_bound
            else
                ps.urow[i] = new_bound
            end
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
    a = abs.(row)
    # building set for heuristicly finding s
    set = build_heuristic_set(a)

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
