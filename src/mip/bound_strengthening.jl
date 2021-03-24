@doc raw"""
    StrengthenSingleIntegerBound <: AbstractRule

Strengthen the bounds of all integer variables using a single constraint.

Let the constraints in the j-th row be
lrow ⩽ rowₛ xₛ + rowⱼ xⱼ ⩽ urow

Denote the U(rowₛ xₛ) be the upper bound of rowₛ xₛ and
L(rowₛ xₛ) be the lower bound of rowₛ xₛ.

Then, lrow - U(rowₛ xₛ) ⩽ rowⱼ xⱼ ⩽ urow - L(rowₛ xₛ).

Thus, if rowⱼ > 0,
xⱼ ⩽ floor((urow - L(rowₛ xₛ)) / rowⱼ) and
xⱼ ⩾ ceil((lrow - U(rowₛ xₛ) xⱼ) / rowⱼ).

If rowⱼ < 0,
xⱼ ⩾ ceil((urow - L(rowₛ xₛ)) / rowⱼ) and
xⱼ ⩽ floor((lrow - U(rowₛ xₛ) xⱼ) / rowⱼ).
"""

struct StrengthenSingleIntegerBound <: AbstractRule
    j::Int
end

function apply!(
    ps::PresolveData{T},
    r::StrengthenSingleIntegerBound,
    config::PresolveOptions{T}
) where {T}
    j = r.j
    ϵ_int = config.IntegerTolerance

    # Column was already removed or the variable is continous.
    ps.colflag[j] || return nothing
    (ps.var_types[j] == CONTINUOUS) && return nothing
    col = ps.pb0.acols[j]
    for (i, row_j) in zip(col.nzind, col.nzval)
        # Row was already removed.
        ps.rowflag[i] || continue
        row = ps.pb0.arows[i]

        lrow = ps.lrow[i]
        urow = ps.urow[i]
        upper = calc_upper_bound_except_one(row, ps.lcol, ps.ucol, j)
        lower = calc_lower_bound_except_one(row, ps.lcol, ps.ucol, j)
        if (row_j > T(0))
            if isfinite(urow) && isfinite(lower)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((urow - lower) / row_j, ϵ_int)))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((lrow - upper) / row_j, ϵ_int)))
            end
        else
            if isfinite(urow) && isfinite(lower)
                ps.lcol[j] = T(max(ps.lcol[j],
                                approx_ceil((urow - lower) / row_j, ϵ_int)))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.ucol[j] = T(min(ps.ucol[j],
                                approx_floor((lrow - upper) / row_j, ϵ_int)))
            end
        end
    end
    return nothing
end


@doc raw"""
    StrengthenSingleContinuousBound <: AbstractRule

Strengthen the bounds of continuous variables using single constraint.

Let the constraints in the j-th row be
lrow ⩽ rowₛ xₛ + rowⱼ xⱼ ⩽ urow

Denote the U(rowₛ xₛ) be the upper bound of rowₛ xₛ and
L(rowₛ xₛ) be the lower bound of rowₛ xₛ.

Then, lrow - U(rowₛ xₛ) ⩽ rowⱼ xⱼ ⩽ urow - L(rowₛ xₛ).

Thus, if rowⱼ > 0,
xⱼ ⩽ (urow - L(rowₛ xₛ)) / rowⱼ and
xⱼ ⩾ (lrow - U(rowₛ xₛ) xⱼ) / rowⱼ.

If rowⱼ < 0,
xⱼ ⩾ (urow - L(rowₛ xₛ)) / rowⱼ and
xⱼ ⩽ (lrow - U(rowₛ xₛ) xⱼ) / rowⱼ.
"""

struct StrengthenSingleContinuousBound <: AbstractRule
    j::Int
end

function apply!(
    ps::PresolveData{T},
    r::StrengthenSingleContinuousBound,
    config::PresolveOptions{T}
) where {T}
    j = r.j

    # Column was already removed or the variable isn't continous.
    ps.colflag[j] || return nothing
    (ps.var_types[j] == CONTINUOUS) || return nothing
    col = ps.pb0.acols[j]
    for (i, row_j) in zip(col.nzind, col.nzval)
        # Row was already removed.
        ps.rowflag[i] || continue
        row = ps.pb0.arows[i]

        lrow = ps.lrow[i]
        urow = ps.urow[i]
        upper = calc_upper_bound_except_one(row, ps.lcol, ps.ucol, j)
        lower = calc_lower_bound_except_one(row, ps.lcol, ps.ucol, j)
        if (row_j > T(0))
            if isfinite(urow) && isfinite(lower)
                ps.ucol[j] = T(min(ps.ucol[j],
                                (urow - lower) / row_j))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.lcol[j] = T(max(ps.lcol[j],
                                (lrow - upper) / row_j))
            end
        else
            if isfinite(urow) && isfinite(lower)
                ps.lcol[j] = T(max(ps.lcol[j],
                                (urow - lower) / row_j))
            end
            if isfinite(lrow) && isfinite(upper)
                ps.ucol[j] = T(min(ps.ucol[j],
                                (lrow - upper) / row_j))
            end
        end
    end
    return nothing
end

"""
    StrengthenBounds <: AbstractRule
Strengthen the bounds of all the integer variables in the problem.

Note: The behavior of `StrengthenBounds` depends on
the ordering of the variables.
"""
struct StrengthenBounds <: AbstractRule end

function apply!(
    ps::PresolveData{T},
    ::StrengthenBounds,
    config::PresolveOptions{T}
) where {T}
    # The problem is LP.
    ps.pb0.is_continuous && return nothing

    for j in 1:ps.pb0.nvar
        if (ps.var_types[j] == CONTINUOUS)
            apply!(ps, StrengthenSingleContinuousBound(j), config)
        else
            apply!(ps, StrengthenSingleIntegerBound(j), config)
        end
    end
    return nothing
end

function calc_upper_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        if (val > 0)
            isfinite(ucol[ind]) || return T(Inf)
            bound += val * ucol[ind]
        else
            isfinite(lcol[ind]) || return T(Inf)
            bound += val * lcol[ind]
        end
    end
    return bound
end

function calc_lower_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        if (val > 0)
            isfinite(lcol[ind]) || return T(-Inf)
            bound += val * lcol[ind]
        else
            isfinite(ucol[ind]) || return T(-Inf)
            bound += val * ucol[ind]
        end
    end
    return bound
end
