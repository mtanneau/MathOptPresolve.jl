"""
Perform Chavatal-Gomory strengthening of inequalities

We apply the Chvatal-Goromory procedure to constraints that
satisfy the following criteria:
1. All variables with nonzero coefficient must be integral
2. The constraints are not ranged
3. If the constraint is of the form ax ⩾ b, then every
   positive coefficient variable must have finite lowerbound and
   negative coefficient variables must have finite upperbound.

   If the constraint is of the from ax ⩽ b, then every
   positive coefficient variable must have finite upperbound and
   negative coefficient variables must have finite lowerbound.

For example, let ax ⩾ b be a constraint where aᵢ ⩾ 0,
and every non-zeoro coefficient in a corresponds to an
integer variable. If there exists an s > 0 such that
⌈aⱼ s⌉ b / ⌈b s⌉ ⩽ aⱼ
and there is k such that
⌈aₖ s⌉ b / ⌈b s⌉ < aₖ (*)
then we will update ax ⩾ b into ⌈a s⌉ x ⩾ ⌈b s⌉

We find such a number s heuristically, i.e we try different
values from the set
{1, t/a_max, t/a_min, (2t-1)/(2a_min)| t = 1,..., r}
where a_max = sup(|aᵢ|), a_min = inf(|aᵢ|) and r is resolution which
is set to be 5 by default

"""
function cg_strengthening!(ps::PresolveData{T}, resolution::Int = 5) where {T}
    ps.pb0.is_continuous && return nothing

    for i in 1:ps.nrow
        ps.rowflag[i] || continue
        nonzero_index = ps.pb0.arows[i].nzind

        if isfinite(ps.urow[i]) && isfinite(ps.lrow[i]) ||
            any(ps.var_types[nonzero_index] .== CONTINUOUS)
            # skip ranged constraints
            # and any constraint that involves continuous variables
            continue
        else
            cg_inequality!(ps, i, resolution)
        end
    end
    col_sync(ps)
    return nothing
end

function cg_inequality!(ps::PresolveData{T}, i::Int, resolution::Int) where {T}
    # perform CG strengthening inequality for a single constraint

    if isfinite(ps.lrow[i]) # var_bound_1 is lower bound
        var_bound_1 = ps.lcol
        var_bound_2 = ps.ucol
        b = ps.lrow[i]
        round_function = ceil
        sign = 1 # indicator for which direction of inequality we're processing
    elseif isfinite(ps.urow[i])
        var_bound_1 = ps.ucol
        var_bound_2 = ps.lcol
        b = ps.urow[i]
        round_function = floor
        sign = -1
    end

    r = ps.pb0.arows[i].nzval
    nonzero_index = ps.pb0.arows[i].nzind
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
    set = build_heuristic_set(a, resolution)

    # perform search for s and update
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

function build_heuristic_set(a::Vector{T}, resolution::Int) where {T}
    # building a heuristic set, for every constraints, we will try every every
    # value in this set to see if any of them satisfies condition (*)
    a_max = maximum(a)
    a_min = minimum(a)

    set = [T(1)]
    for t in 1:resolution
        append!(set, [t/a_max, t/a_min, (2*t-1)/(2*a_min)])
    end
    return unique(set)
end
