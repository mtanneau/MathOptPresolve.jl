@doc raw"""
    FixSingleVariable <: AbstractRule

Eliminate a single variable ``x_{j}`` if fixed.

## Presolve

A variable ``x_{j}`` with lower and upper bounds ``l_{j}, u_{j}``
is fixed to its lower bound if
```math
| l_{j} - u_{j}| \leq ϵ,
```
where ``\epsilon`` is a prescribed tolerance.

The variable is then eliminated from the problem as follows.
Forall rows $i$,
```math
l^{r}_{i} \leq a_{i, j}x_{j} + \sum_{k \neq j} a_{i, k} x_{k} \leq u^{r}_{i}
```
is transformed into
```math
l^{r}_{i} - a_{i, j} l_{j} \leq \sum_{k \neq j} a_{i, k} x_{k} \leq u^{r}_{i} - a_{i, j} l_{j}.
```
The objective constant is updated to $c_{0} + c_{j} l_{j}$.

### Integer infeasibility

The problem is infeasible if an integer variable is fixed to a fractional value, namely,
if the fixed variable ``x_{j}`` is integer (or implied integer) and
```math
\min(l_{j} - ⌊ l_{j} ⌋, ⌈ l_{j} ⌉ - l_{j}) \geq \epsilon_{int},
```
where ``\epsilon_{int}`` is a prescribed tolerance.

If ``x_{j}`` is binary, then we also check whether ``l_{j} \in \{0, 1\}`` --
up to the prescribed tolerance ``\epsilon_{int}``, otherwise the problem is infeasible.

## Postsolve

For dual variables, first compute
```math
z_{j} = c_{j} - \sum_{i} a_{i, j} (y_{i}^{l} - y_{i}^{u}),
```
then recover ``z_{j}^{l} = z_{j}^{+}`` and ``z_{j}^{u} = z_{j}^{-}``.


## Misc

* This is a primal reduction.
* This reduction is non-destructive.
* This reduction does not create any fill-in.
"""
struct FixSingleVariable <: AbstractRule
    j::Int
end

@doc raw"""
    FixedVariable{T} <: AbstractReduction{T}

Variable ``x_{j}`` has been fixed

See [`FixSingleVariable`](@ref)

"""
struct FixedVariable{T} <: AbstractReduction{T}
    j::Int  # variable index
    x::T  # primal value
    c::T  # current objective coeff
    col::Col{T}  # current column
end

function apply!(
    ps::PresolveData{T},
    r::FixSingleVariable,
    config::PresolveOptions{T}
) where {T}
    j = r.j
    ϵ = config.PrimalTolerance
    ϵ_int = config.IntegerTolerance

    ps.colflag[j] || return nothing  # Column was already removed

    # Get current bounds
    lb, ub = ps.lcol[j], ps.ucol[j]

    abs(ub - lb) <= ϵ || return nothing

    @debug "Fixing variable $j to $lb"

    # Infeasibility check for integer variables
    vt = ps.pb0.var_types[j]
    if is_integer(vt)
        # Check if fixing to fractional value
        f = min(lb - floor(lb), ceil(lb) - lb)
        if f >= ϵ_int || (vt == BINARY
            && ((lb <= -ϵ_int) || (lb >= 1 + ϵ_int))
        )
            # MIP problem is infeasible
            ps.status = PRIMAL_INFEASIBLE
            return nothing
        end
    end

    col = ps.pb0.acols[j]
    cj = ps.obj[j]

    # Remove column
    ps.colflag[j] = false
    ps.ncol -= 1
    ps.updated = true

    # TODO: make this more efficient
    push!(ps.ops, FixedVariable(j, lb, cj, Col(
        [i for i in col.nzind if ps.rowflag[i]],
        [aij for (i, aij) in zip(col.nzind, col.nzval) if ps.rowflag[i]]
    )))

    # Update objective constant
    ps.obj0 += cj * lb

    # Update rows
    # TODO: update list of row singletons and empty rows
    for (i, aij) in zip(col.nzind, col.nzval)
        ps.rowflag[i] || continue  # This row is no longer in the problem
        iszero(aij) && continue  # Skip if coefficient is zero

        # Update row bounds
        ps.lrow[i] -= aij * lb
        ps.urow[i] -= aij * lb

        ps.nzrow[i] -= 1

        # Check row
        if ps.nzrow[i] == 0
            remove_empty_row!(ps, i)
        elseif ps.nzrow == 1
            push!(ps.row_singletons, i)
        end
    end  # row update

    # Done
    return nothing
end

function postsolve!(sol::Solution{T}, r::FixedVariable{T}) where{T}
    sol.x[r.j] = r.x
    s = sol.is_dual_ray ? zero(T) : r.c
    for (i, aij) in zip(r.col.nzind, r.col.nzval)
        s -= aij * (sol.y_lower[i] - sol.y_upper[i])
    end
    sol.s_lower[r.j] = pos_part(s)
    sol.s_upper[r.j] = neg_part(s)
    return nothing
end

"""
    RemoveFixedVariables <: AbstractRule

Remove all fixed variables in the problem.

See [`FixSingleVariable`](@ref).
"""
struct RemoveFixedVariables <: AbstractRule end

function apply!(
    ps::PresolveData{T},
    ::RemoveFixedVariables,
    config::PresolveOptions{T}
) where {T}

    for (j, flag) in enumerate(ps.colflag)
        apply!(ps, FixSingleVariable(j), config)
        ps.status == NOT_INFERRED || return nothing
    end

    return nothing
end
