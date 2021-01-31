@doc raw"""
    RemoveEmptyColumn <: AbstractRule

Eliminate a single column if it is empty, i.e. all coefficients are zero.

## Presolve

If column $j$ is empty, variable $x_{j}$ is fixed to $x_{j}^{*}$ as follows

| $c_{j}$ |        $l_{j}$ |        $u_{j}$ | $x^{*}_{j}$ |
|--------:|---------------:|---------------:|------------:|
|    $<0$ |              - | $∈ \mathbb{R}$ |     $u_{j}$ |
|    $<0$ |              - |           $+∞$ |    $+∞^{†}$ |
|     $0$ | $∈ \mathbb{R}$ |              - |     $l_{j}$ |
|     $0$ | $-∞$           | $∈ \mathbb{R}$ |     $u_{j}$ |
|     $0$ | $-∞$           |           $+∞$ |     $0$     |
|    $>0$ | $∈ \mathbb{R}$ |              - |     $l_{j}$ |
|    $>0$ | $-∞$           |              - |    $-∞^{†}$ |
``^{†}`` problem is dual infeasible

## Dual infeasibility (unboundedness) checks

For a prescribed tolerance $ϵ ≥ 0$, the problem is dual infeasible if
    column $j$ is empty and either
1. Variable $j$ has no lower bound, and its objective coefficient is positive,
    i.e., $l^{c}_{j} = -∞$ and $c_{j} > ϵ$.
2. column $j$ is empty, has no upper bound, and its objective coefficient is negative,
    i.e., $l^{c}_{j} = +∞$ and $c_{j} < -ϵ$.

The corresponding unbounded rays (dual infeasibility certificates) are
1.  $x = -e_{j}$,
2.  $x = e_{j}$,

where $e_{j}$ is the $j$-th vector of the canonical basis.

## Postsolve

TODO

## Misc

* This is a dual reduction, and it can handle integer variables.
* This reduction is non-destructive.
* This reduction does not create any fill-in.
"""
struct RemoveEmptyColumn <: AbstractRule
    j::Int
end

@doc raw"""
    EmptyColumn{T} <: AbstractReduction{T}

Column $j$ is empty (all coefficients are zero).
"""
struct EmptyColumn{T} <: AbstractReduction{T}
    j::Int  # variable index
    x::T  # Variable value
    s::T  # Reduced cost
end

function apply!(
    ps::PresolveData{T},
    r::RemoveEmptyColumn,
    config::PresolveOptions{T},
) where {T}
    j = r.j
    ϵ = config.DualTolerance

    ps.pb0.is_continuous || error("Empty column routine currently only supported for LPs.")

    # Sanity check
    (ps.colflag[j] && (ps.nzcol[j] == 0)) || return nothing

    # Remove column
    lb, ub = ps.lcol[j], ps.ucol[j]
    cj = ps.obj[j]
    @debug "Removing empty column $j" cj lb ub

    if cj > zero(T)
        if isfinite(lb)
            # Set variable to lower bound
            # Update objective constant
            ps.obj0 += lb * cj
            push!(ps.ops, EmptyColumn(j, lb, cj))
        else
            # Problem is dual infeasible
            @debug "Column $j is (lower) unbounded"
            ps.status = DUAL_INFEASIBLE
            ps.updated = true

            # Resize problem
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(T)
            ps.solution.y_lower .= zero(T)
            ps.solution.y_upper .= zero(T)
            ps.solution.s_lower .= zero(T)
            ps.solution.s_upper .= zero(T)

            # Unbounded ray: xj = -1
            ps.solution.primal_status = INFEASIBILITY_CERTIFICATE
            ps.solution.dual_status = NO_SOLUTION
            ps.solution.is_primal_ray = true
            ps.solution.is_dual_ray = false
            ps.solution.z_primal = ps.solution.z_dual = -T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.x[j_] = -one(T)

            return nothing
        end
    elseif cj < zero(T)
        if isfinite(ub)
            # Set variable to upper bound
            # Update objective constant
            ps.obj0 += ub * cj
            push!(ps.ops, EmptyColumn(j, ub, cj))
        else
            # Problem is dual infeasible
            @debug "Column $j is (upper) unbounded"
            ps.status = DUAL_INFEASIBLE
            ps.updated = true

            # Resize problem
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(T)
            ps.solution.y_lower .= zero(T)
            ps.solution.y_upper .= zero(T)
            ps.solution.s_lower .= zero(T)
            ps.solution.s_upper .= zero(T)

            # Unbounded ray: xj = 1
            ps.solution.primal_status = INFEASIBILITY_CERTIFICATE
            ps.solution.dual_status = NO_SOLUTION
            ps.solution.is_primal_ray = true
            ps.solution.is_dual_ray = false
            ps.solution.z_primal = ps.solution.z_dual = -T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.x[j_] = one(T)

            return
        end
    else
        # Any feasible value works
        if isfinite(lb)
            push!(ps.ops, EmptyColumn(j, lb, zero(T)))
        elseif isfinite(ub)
            push!(ps.ops, EmptyColumn(j, ub, zero(T)))
        else
            # Free variable with zero coefficient and empty column
            push!(ps.ops, EmptyColumn(j, zero(T), zero(T)))
        end

    end

    # Book=keeping
    ps.colflag[j] = false
    ps.updated = true
    ps.ncol -= 1
    return nothing
end

function postsolve!(sol::Solution{T}, op::EmptyColumn{T}) where {T}
    sol.x[op.j] = op.x
    sol.s_lower[op.j] = pos_part(op.s)
    sol.s_upper[op.j] = neg_part(op.s)
    return nothing
end

"""
    RemoveEmptyColumns <: AbstractRule

Remove all empty columns in the problem.

See [`RemoveEmptyColumn`](@ref).
"""
struct RemoveEmptyColumns <: AbstractRule end

function apply!(
    ps::PresolveData{T},
    ::RemoveEmptyColumns,
    config::PresolveOptions{T}
) where {T}
    for j in 1:ps.pb0.nvar
        if ps.colflag[j] && (ps.nzcol[j] == 0)
            apply!(ps, RemoveEmptyColumn(j), config)
            ps.status == NOT_INFERRED || return nothing
        end
    end
    return nothing
end
