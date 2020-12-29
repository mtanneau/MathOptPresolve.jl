@doc raw"""
    RemoveEmptyRow <: AbstractRule

Eliminate a single row if it is empty, i.e. all coefficients are zero.

In particular, consider row $i$ given as
```math
l^{r}_{i} ≤ \sum_{j} a_{i, j} x_{j} ≤ u^{r}_{i}.
```
We eliminate row $i$ if $a_{i,j} = 0$ for each variable index $j$.

## Presolve

If row $i$ is empty, the corresponding constraint is completely removed from
the problem.

## Infeasibility checks

For a prescribed tolerance $ϵ ≥ 0$, the problem is infeasible if either
1. row $i$ is empty and its lower bound is positive, i.e., $l^{r}_{i} > ϵ$, or
2. row $i$ is empty and its upper bound is negative, i.e., $u^{r}_{i} < -ϵ$.

The corresponding Farkas infeasibilities certificates are
1.  $y^{l} = e_{i}$, $y^{u} = 0$, $z^{l} = z^{u} = 0$,
2.  $y^{u} = e_{i}$, $y^{l} = 0$, $z^{l} = z^{u} = 0$,

where $e_{i}$ is the $i$-th vector of the canonical basis.

## Postsolve

If the problem is feasible, $y^{l}_{i} = y^{u}_{i} = 0$.

If the problem is infeasible, see the Farkas rays above.

## Misc

* This is a primal reduction.
* This reduction is non-destructive.
* This reduction does not create any fill-in.
"""
struct RemoveEmptyRow <: AbstractRule
    i::Int
end

@doc raw"""
    EmptyRow{T} <: AbstractReduction{T}

Row $i$ is empty (all variable coefficients are zero).
"""
struct EmptyRow{T} <: AbstractReduction{T}
    i::Int  # row index
    y::T  # dual multiplier
end

function apply!(
    ps::PresolveData{T},
    r::RemoveEmptyRow,
    config::PresolveOptions{T}
) where {T}
    i = r.i
    ϵ = config.PrimalTolerance

    # If row isn't empty, we can't remove it
    (ps.rowflag[i] && ps.nzrow[i] == 0) || return nothing

    # Check bounds
    lb, ub = ps.lrow[i], ps.urow[i]

    if ub < -ϵ
        # Infeasible
        @debug "Row $i is primal infeasible"
        ps.status = PRIMAL_INFEASIBLE
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(T)
        ps.solution.y_lower .= zero(T)
        ps.solution.y_upper .= zero(T)
        ps.solution.s_lower .= zero(T)
        ps.solution.s_upper .= zero(T)

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = NO_SOLUTION
        ps.solution.dual_status = INFEASIBILITY_CERTIFICATE
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = T(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_upper[i] = one(T)
        return
    elseif lb > ϵ
        @debug "Row $i is primal infeasible"
        ps.status = PRIMAL_INFEASIBLE
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(T)
        ps.solution.y_lower .= zero(T)
        ps.solution.y_upper .= zero(T)
        ps.solution.s_lower .= zero(T)
        ps.solution.s_upper .= zero(T)

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = NO_SOLUTION
        ps.solution.dual_status = INFEASIBILITY_CERTIFICATE
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = T(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_lower[i_] = one(T)
        return
    else
        push!(ps.ops, EmptyRow(i, zero(T)))
    end

    # Book-keeping
    ps.updated = true
    ps.rowflag[i] = false
    ps.nrow -= 1

    return nothing
end

function postsolve!(sol::Solution{T}, r::EmptyRow{T}) where {T}
    sol.y_lower[r.i] = pos_part(r.y)
    sol.y_upper[r.i] = neg_part(r.y)
    return nothing
end

"""
    RemoveEmptyRows <: AbstractRule

Remove all empty rows in the problem.

See [`RemoveEmptyRow`](@ref).
"""
struct RemoveEmptyRows <: AbstractRule end

function apply!(
    ps::PresolveData{T},
    ::RemoveEmptyRows,
    config::PresolveOptions{T}
) where {T}
    for i in 1:ps.pb0.ncon
        if ps.rowflag[i] && (ps.nzrow[i] == 0)
            apply!(ps, RemoveEmptyRow(i), config)
            ps.status == NOT_INFERRED || return nothing
        end
    end
    return nothing
end
