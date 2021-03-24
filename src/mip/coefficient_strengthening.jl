"""
Perform coefficient strengthening on each row and on each integer variable.

In particular, we consider the i-th constraints of the form
aᵢ x ⩽ bᵢ and where xⱼ is a integer/binary variable.
If aᵢⱼ ⩾ d = bᵢ - Mᵢⱼ - aᵢⱼ (uⱼ - 1) > 0
where Mᵢⱼ is maximal activity of the i-th row without considering xⱼ
and uⱼ is upper bound of xⱼ, then we transform
aᵢⱼ <- aᵢⱼ - d
bᵢ <- bᵢ - duⱼ .

A similar update rule is applied for the constraints of the form cᵢ ⩽ aᵢ x.

In case the constraint is ranged, i.e, cᵢ ⩽ aᵢ x ⩽ bᵢ, no coefficient strengthening is performed.

The new coefficients are stored by updating ps.pb0 .
If some coefficient is reduced to 0, it is kept as 0 in both ps.pb0.arows and ps.pb0.acols,
and the number of nonzeros in each row/column is updated.

If a row/column is reduced to 0 after this step, it will be set to be inactive.
"""

function upperbound_strengthening(ps::PresolveData{T}, i::Int, j::Int, coef, max_act::T) where {T}
    # perform coef strengthening for one constraint of the form a'x <= u
    new_bound = ps.urow[i]
    new_coef = coef
    if coef > 0
        d = new_bound - max_act - coef * (ps.ucol[j]-1)
        if coef >= d > 0
            new_coef = new_coef - d
            new_bound -= d*ps.ucol[j]
        end
    elseif coef < 0
        d = new_bound - max_act - coef * (ps.lcol[j]+1)
        if -coef >= d > 0
            new_coef = new_coef + d
            new_bound += d*ps.lcol[j]
        end
    end
    return new_coef, new_bound
end

function lowerbound_strengthening(ps::PresolveData{T}, i::Int, j::Int, coef, min_act::T) where {T}
    # perform coef strengthening for one constraint of the form l < = a'x
    new_bound = ps.lrow[i]
    new_coef = coef
    if coef > 0
        d = -new_bound + min_act + coef*(ps.lcol[j]+1)
        if coef >= d > 0
            new_coef = coef - d
            new_bound -= d*ps.lcol[j]
        end
    elseif coef < 0
        d = -new_bound + min_act + coef*(ps.ucol[j]-1)
        if -coef >= d > 0
            new_coef = coef + d
            new_bound += d*ps.ucol[j]
        end
    end
    return new_coef, new_bound
end

"""
    coefficient_strengthening!(ps::PresolveData)

Perform coefficient strengthening for integer/binary variables
in every constraint.

Called once at the beginning of the presolve procedure.
"""

function coefficient_strengthening!(ps::PresolveData{T}) where {T}
    # perform coefficient stregthening but if there is a coefficient that is reduced to 0
    # it is explicitly stored in the ps.pb0.arows and ps.pb0.acols

    ps.pb0.is_continuous && return nothing

    # keep track of index for each var to update ps.acols
    # use this to find which index of ps.pb0.acols[i].nzval to update
    i_index = zeros(Int, ps.pb0.nvar)

    for i in 1:ps.pb0.ncon
        ps.rowflag[i] || continue

        lrow = ps.lrow[i]
        urow = ps.urow[i]
        if isfinite(lrow) && isfinite(urow) #skip ranged constraints
            continue
        end

        row = ps.pb0.arows[i]
        if all(ps.var_types[row.nzind] .== CONTINUOUS) # at least 1 integer
            continue
        end

        sup = maximal_activity(ps, i)
        inf = minimal_activity(ps, i)

        j_index = 0 # keep track of index of variable j in row.nzind & row.nzval
        for j in row.nzind
            j_index += 1
            i_index[j] += 1
            if ps.var_types[j] == CONTINUOUS || !ps.colflag[j]
                continue
            else
                coef = row.nzval[j_index]

                if isfinite(urow)
                    if coef > 0
                        max_act = sup - coef * ps.ucol[j] # maximal activity of every variables except j
                    else
                        max_act = sup - coef * ps.lcol[j]
                    end
                    new_coef, new_bound = upperbound_strengthening(ps, i, j, coef, max_act)
                    # update problem
                    row.nzval[j_index] = new_coef
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    ps.urow[i] = new_bound
                    #update nonzero
                    if !iszero(coef) && iszero(new_coef)
                        ps.nzrow[i] -= 1
                        ps.nzcol[j] -= 1
                    end
                    #update sup
                    if coef > 0
                        sup -= (coef - new_coef) * ps.ucol[j]
                    else
                        sup -= (coef - new_coef) * ps.lcol[j]
                    end
                elseif isfinite(lrow)
                    if coef > 0
                        min_act = inf - coef * ps.lcol[j] # minimal activity of every variables except j
                    else
                        min_act = inf - coef * ps.ucol[j]
                    end
                    new_coef, new_bound = lowerbound_strengthening(ps, i, j, coef, min_act)
                    #update problem
                    row.nzval[j_index] = new_coef
                    ps.pb0.acols[j].nzval[i_index[j]] = new_coef
                    ps.lrow[i] = new_bound
                    #update nonzero
                    if !iszero(coef) && iszero(new_coef)
                        ps.nzrow[i] -= 1
                        ps.nzcol[j] -= 1
                    end
                    #update inf
                    if coef > 0
                        inf -= (coef - new_coef) * ps.lcol[j]
                    else
                        inf -= (coef - new_coef) * ps.ucol[j]
                    end
                end
            end
        end
    end
    return nothing
end
