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

For example, let ax ⩾ b be a constraint where a ⩾ 0,
and every non-zero coefficient in a corresponds to an
integer variable. If there exists an s > 0 such that
⌈aⱼ s⌉ b / ⌈b s⌉ ⩽ aⱼ
for all j and there is k such that
⌈aₖ s⌉ b / ⌈b s⌉ < aₖ (*)
then we will update ax ⩾ b into ⌈a s⌉ x ⩾ ⌈b s⌉.

We find such a number s heuristically, i.e we try different
values from the set
{1, t/a_max, t/a_min, (2t-1)/(2a_min)| t = 1,..., r}
where a_max = max(|aᵢ|), a_min = min(|aᵢ|) and r is resolution which
is set to be 5 by default.

When performing CG strengthening, we round up positive
coefficients and round down negative ones, so this procedure
will not reduce any coefficient to 0.
"""

function cg_strengthening!(ps::PresolveData{T}, resolution::Int = 5) where {T}
    ps.pb0.is_continuous && return nothing
    for i in 1:ps.pb0.ncon
        ps.rowflag[i] || continue
        nzind = ps.pb0.arows[i].nzind
        nzval = ps.pb0.arows[i].nzval
        finite_lb = isfinite(ps.lrow[i])
        finite_ub = isfinite(ps.urow[i])

        if finite_lb && finite_ub ||
            any(ps.var_types[[j for j in nzind if ps.colflag[j]]] .== CONTINUOUS)
            # skip ranged constraints
            # and any constraint that involves continuous variables
            continue
        elseif finite_lb
            s, newbound = lowerbound_cg_inequality(nzval, nzind, ps.lrow[i],
                                    ps.colflag, ps.lcol, ps.ucol, resolution)
            if s != 0
                ps.lrow[i] = newbound
                ps.pb0.arows[i].nzval = [aij > 0 ? ceil(s*aij) : floor(s*aij) for aij in nzval]
            end
        elseif finite_ub
            s, newbound = lowerbound_cg_inequality(-nzval, nzind, -ps.urow[i],
                                    ps.colflag, ps.lcol, ps.ucol, resolution)
            if s != 0
                ps.urow[i] = -newbound
                ps.pb0.arows[i].nzval = [aij > 0 ? ceil(s*aij) : floor(s*aij) for aij in nzval]
            end
        else
            error("Unbounded Constraint!")
        end
    end
    return nothing
end

"""
Find the CG strengthening inequality for constraints of the form ax ⩾ b

We return the value s for updating a and the new lower bound. If s is 0, then
no CG inequality is found with the current resolution.
"""
function lowerbound_cg_inequality(
    nzval::Vector{T},
    nzind::Vector{Int},
    b::T,
    colflag::Vector{Bool},
    var_lb::Vector{T},
    var_ub::Vector{T},
    resolution::Int) where {T}

    new_bound = b
    row_pos_index = Int[]
    row_neg_index = Int[]
    col_pos_index = Int[]
    col_neg_index = Int[]
    for (k, (j, aij)) in enumerate(zip(nzind, nzval))
        if aij > 0 && isfinite(var_lb[j])
            if colflag[j]
                new_bound -= aij * var_lb[j]
                push!(row_pos_index, k)
                push!(col_pos_index, j)
            end
        elseif aij < 0 && isfinite(var_ub[j])
            if colflag[j]
                new_bound -= aij * var_ub[j]
                push!(row_neg_index, k)
                push!(col_neg_index, j)
            end
        else
            # if any positive variable has -inf lower bound
            # or any negative variable has inf upper bound, then terminate.
            return 0, new_bound
        end
    end

    # building set for heuristicly finding s.
    a = [abs(aij) for (j, aij) in zip(nzind, nzval) if colflag[j]]
    set = build_heuristic_set(a, resolution)

    # perform search for s and update.
    s = 0
    for h in set
        found_s = false
        for val in a
            if val*ceil(new_bound*h) < ceil(h*val)*new_bound
                found_s = false
                break
            elseif val*ceil(new_bound*h) > ceil(h*val)*new_bound
                found_s = true
            end
        end

        if found_s
            s = h
            new_bound = ceil(s*new_bound) +
            dot(ceil.(s*a[row_pos_index]), var_lb[col_pos_index]) -
            dot(ceil.(s*a[row_neg_index]), var_ub[col_neg_index])
            break
        end
    end
    return s, new_bound
end


function build_heuristic_set(a::Vector{T}, resolution::Int) where {T}
    # building a heuristic set, for every constraints, we will try every
    # value in this set to see if any of them satisfies condition (*).
    a_min, a_max = extrema(a)
    set = [T(1)]
    for t in 1:resolution
        append!(set, [t/a_max, t/a_min, (2*t-1)/(2*a_min)])
    end
    return unique(set)
end
