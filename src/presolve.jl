"""
    AbstractReduction{T}

Abstract type for presolve transformations.
"""
abstract type AbstractReduction{T} end

abstract type AbstractRule end

"""
    apply!(ps, rule, config) -> Nothing

Apply rule `rule`.

# Arguments
* `ps::PresolveData{T}`
* `rule::AbstractRule`
* `config::PresolveOptions{T}`
"""
function apply! end

"""
    PresolveData{T}

Stores information about an LP in the form
```
    min     cᵀx + c₀
    s.t.    lr ⩽ Ax ⩽ ur
            lc ⩽  x ⩽ uc
```
whose dual writes
```
    max     lrᵀy⁺ - urᵀy⁻ + lcᵀs⁺ - ucᵀs⁻
    s.t.     Aᵀy⁺ -  Aᵀy⁻ +    s⁺ -    s⁻ = c
               y⁺,     y⁻,     s⁺,     s⁻ ⩾ 0
```
"""
mutable struct PresolveData{T}
    updated::Bool
    status::ModelStatus

    # Original problem
    pb0::ProblemData{T}
    # Reduced problem
    # Nothing until the reduced problem is extracted
    pb_red::Union{Nothing,ProblemData{T}}
    solution::Solution{T}  # only used if presolve solves the problem

    # Presolved data

    # Active rows and columns
    rowflag::Vector{Bool}
    colflag::Vector{Bool}

    # Non-zeros in rows and columns
    nzrow::Vector{Int}
    nzcol::Vector{Int}

    # Objective
    objsense::Bool
    obj::Vector{T}
    obj0::T

    # Current number of constraints/variables in presolved problem
    nrow::Int
    ncol::Int

    # Primal bounds
    lrow::Vector{T}
    urow::Vector{T}
    lcol::Vector{T}
    ucol::Vector{T}

    # Dual bounds
    ly::Vector{T}
    uy::Vector{T}
    ls::Vector{T}
    us::Vector{T}

    # Variable types
    var_types::Vector{VariableType}

    # Scaling
    row_scaling::Vector{T}
    col_scaling::Vector{T}
    # TODO: objective and RHS scaling

    # Old <-> new index mapping
    # Instantiated only after pre-solve is performed
    new_con_idx::Vector{Int}
    new_var_idx::Vector{Int}
    old_con_idx::Vector{Int}
    old_var_idx::Vector{Int}

    # Singletons
    row_singletons::Vector{Int}  # Row singletons
    free_col_singletons::Vector{Int}  # (implied) free column singletons

    # TODO: set of transformations for pre-post crush
    ops::Vector{AbstractReduction{T}}

    function PresolveData(pb::ProblemData{T}) where {T}
        ps = new{T}()

        ps.updated = false
        ps.status = NOT_INFERRED

        ps.pb0 = pb
        ps.pb_red = nothing
        ps.solution = Solution{T}(pb.ncon, pb.nvar)

        ps.nrow = pb.ncon
        ps.ncol = pb.nvar

        # All rows and columns are active
        ps.rowflag = trues(ps.nrow)
        ps.colflag = trues(ps.ncol)

        # Number of non-zeros in rows/columns
        ps.nzrow = zeros(Int, ps.nrow)
        ps.nzcol = zeros(Int, ps.ncol)
        for (j, col) in enumerate(pb.acols)
            for (i, aij) in zip(col.nzind, col.nzval)
                ps.nzcol[j] += !iszero(aij)
                ps.nzrow[i] += !iszero(aij)
            end
        end

        # Objective
        ps.objsense = pb.objsense
        if pb.objsense
            ps.obj  = copy(pb.obj)
            ps.obj0 = pb.obj0
        else
            # Maximization problem: negate the objective for pre-solve
            # This will be undone when extracting the reduced problem
            ps.obj  = -copy(pb.obj)
            ps.obj0 = -pb.obj0
        end

        # Copy primal bounds
        ps.lrow = copy(pb.lcon)
        ps.urow = copy(pb.ucon)
        ps.lcol = copy(pb.lvar)
        ps.ucol = copy(pb.uvar)

        # Set dual bounds
        ps.ly = Vector{T}(undef, ps.nrow)
        ps.uy = Vector{T}(undef, ps.nrow)
        ps.ls = Vector{T}(undef, ps.ncol)
        ps.us = Vector{T}(undef, ps.ncol)
        for (i, (lc, uc)) in enumerate(zip(ps.lrow, ps.urow))
            ps.ly[i] = (uc == T(Inf)) ? zero(T) : T(-Inf)
            ps.uy[i] = (lc == T(-Inf)) ? zero(T) : T(Inf)
        end
        for (j, (lv, uv)) in enumerate(zip(ps.lcol, ps.ucol))
            ps.ls[j] = (uv == T(Inf)) ? zero(T) : T(-Inf)
            ps.us[j] = (lv == T(-Inf)) ? zero(T) : T(Inf)
        end

        # Variable types
        ps.var_types = copy(pb.var_types)

        # Scalings
        ps.row_scaling = ones(T, ps.nrow)
        ps.col_scaling = ones(T, ps.ncol)

        # Index mappings
        ps.new_con_idx = Int[]
        ps.new_var_idx = Int[]
        ps.old_con_idx = Int[]
        ps.old_var_idx = Int[]

        # Singletons
        ps.row_singletons = Int[]
        ps.free_col_singletons = Int[]

        ps.ops = AbstractReduction{T}[]

        return ps
    end
end

# Extract pre-solved problem data, to be passed to the IPM solver
function extract_reduced_problem!(ps::PresolveData{T}) where {T}

    pb = ProblemData{T}()

    pb.ncon = sum(ps.rowflag)
    pb.nvar = sum(ps.colflag)

    pb.objsense = ps.objsense
    if pb.objsense
        pb.obj0 = ps.obj0
        pb.obj = ps.obj[ps.colflag]
    else
        pb.obj0 = -ps.obj0
        pb.obj = -ps.obj[ps.colflag]
    end

    # Primal bounds
    pb.lvar = ps.lcol[ps.colflag]
    pb.uvar = ps.ucol[ps.colflag]
    pb.lcon = ps.lrow[ps.rowflag]
    pb.ucon = ps.urow[ps.rowflag]

    # Extract new rows
    pb.arows = Vector{Row{T}}(undef, pb.ncon)
    inew = 0
    for (iold, row) in enumerate(ps.pb0.arows)
        ps.rowflag[iold] || continue

        inew += 1
        # Compute new row
        rind = Vector{Int}(undef, ps.nzrow[iold])
        rval = Vector{T}(undef, ps.nzrow[iold])

        k = 0
        for (jold, aij) in zip(row.nzind, row.nzval)
            ps.colflag[jold] || continue
            iszero(aij) && continue

            # Set new coefficient
            k += 1
            rind[k] = ps.new_var_idx[jold]
            rval[k] = aij
        end

        # Set new row
        pb.arows[inew] = Row{T}(rind, rval)
    end

    # Extract new columns
    pb.acols = Vector{Col{T}}(undef, pb.nvar)
    pb.var_types = Vector{VariableType}(undef, pb.nvar)
    jnew = 0
    for (jold, col) in enumerate(ps.pb0.acols)
        ps.colflag[jold] || continue

        jnew += 1
        # Compute new row
        cind = Vector{Int}(undef, ps.nzcol[jold])
        cval = Vector{T}(undef, ps.nzcol[jold])

        k = 0
        for (iold, aij) in zip(col.nzind, col.nzval)
            ps.rowflag[iold] || continue
            iszero(aij) && continue

            # Set new coefficient
            k += 1
            cind[k] = ps.new_con_idx[iold]
            cval[k] = aij
        end

        # Set new column
        pb.acols[jnew] = Col{T}(cind, cval)

        pb.var_types[jnew] = ps.var_types[jold]
    end

    # Scaling
    rscale = zeros(T, ps.nrow)
    cscale = zeros(T, ps.ncol)

    # Compute norm of each row and column
    # TODO: use a parameter p and do norm(.., p)
    p = 2
    for (i, row) in enumerate(pb.arows)
        r = norm(row.nzval, p)
        rscale[i] = r > zero(T) ? r : one(T)
    end
    for (j, col) in enumerate(pb.acols)
        r = norm(col.nzval, p)
        cscale[j] = r > zero(T) ? r : one(T)
    end

    map!(sqrt, cscale, cscale)
    map!(sqrt, rscale, rscale)

    # Rows
    for (i, row) in enumerate(pb.arows)
        # Scale row coefficients
        for (k, j) in enumerate(row.nzind)
            row.nzval[k] /= (rscale[i] * cscale[j])
        end
        # Scale row bounds
        pb.lcon[i] /= rscale[i]
        pb.ucon[i] /= rscale[i]
    end
    # Columns
    for (j, col) in enumerate(pb.acols)
        # Scale column coefficients
        for (k, i) in enumerate(col.nzind)
            col.nzval[k] /= (rscale[i] * cscale[j])
        end
        # Scale objective and variable bounds
        pb.obj[j]  /= cscale[j]
        pb.lvar[j] *= cscale[j]
        pb.uvar[j] *= cscale[j]
    end

    # Record scaling
    @debug "Scaling info" extrema(rscale) extrema(cscale)
    ps.row_scaling = rscale
    ps.col_scaling = cscale

    # Done
    ps.pb_red = pb
    return nothing
end

include("lp/empty_row.jl")
include("lp/empty_column.jl")
include("lp/fixed_variable.jl")
include("lp/row_singleton.jl")
include("lp/forcing_row.jl")
include("lp/free_column_singleton.jl")
include("lp/dominated_column.jl")

include("mip/round_integer_bounds.jl")
include("mip/bound_strengthening.jl")


"""
    postsolve!

Perform post-solve.
"""
function postsolve!(sol::Solution{T}, sol_::Solution{T}, ps::PresolveData{T}) where {T}

    # Check dimensions
    (sol_.m, sol_.n) == (ps.nrow, ps.ncol) || error(
        "Inner solution has size $((sol_.m, sol_.n)) but presolved problem has size $((ps.nrow, ps.ncol))"
    )
    (sol.m, sol.n) == (ps.pb0.ncon, ps.pb0.nvar) || error(
        "Solution has size $((sol.m, sol.n)) but original problem has size $((ps.pb0.ncon, ps.pb0.nvar))"
    )

    # Copy solution status and objective values
    sol.primal_status = sol_.primal_status
    sol.dual_status = sol_.dual_status
    sol.is_primal_ray = sol_.is_primal_ray
    sol.is_dual_ray = sol_.is_dual_ray
    sol.z_primal = sol_.z_primal
    sol.z_dual = sol_.z_dual

    # Extract and un-scale inner solution components
    # TODO: create a AbstractReduction for scaling
    for (j_, j) in enumerate(ps.old_var_idx)
        sol.x[j] = sol_.x[j_] / ps.col_scaling[j_]
        sol.s_lower[j] = sol_.s_lower[j_] * ps.col_scaling[j_]
        sol.s_upper[j] = sol_.s_upper[j_] * ps.col_scaling[j_]
    end
    for (i_, i) in enumerate(ps.old_con_idx)
        sol.y_lower[i] = sol_.y_lower[i_] / ps.row_scaling[i_]
        sol.y_upper[i] = sol_.y_upper[i_] / ps.row_scaling[i_]
    end

    # Reverse transformations
    for op in Iterators.reverse(ps.ops)
        postsolve!(sol, op)
    end

    # Compute row primals
    for (i, row) in enumerate(ps.pb0.arows)
        sol.Ax[i] = zero(T)
        for (j, aij) in zip(row.nzind, row.nzval)
            sol.Ax[i] += aij * sol.x[j]
        end
    end

    # Done
    return nothing
end

"""
Helper macro for running a presolve routine and returning if
that routine was able to infer the model status (e.g. optimal).
"""
macro _return_if_inferred(expr)
    @assert expr.head == :call
    @assert length(expr.args) >= 2
    if length(expr.args) == 2
        func = expr.args[1]
        ps = expr.args[2]
        return esc(quote
            $func($ps)
            $ps.status == NOT_INFERRED || return $ps.status
        end)
    elseif length(expr.args) == 4
        func = expr.args[1]
        ps = expr.args[2]
        rule = expr.args[3]
        config = expr.args[4]
        return esc(quote
            $func($ps, $rule, $config)
            $ps.status == NOT_INFERRED || return $ps.status
        end)
    else
        error("_return_if_inferred needs 2 or 4 arguments")
    end
end

"""
    presolve(pb::ProblemData)

Perform pre-solve.
"""
function presolve!(ps::PresolveData{T}) where {T}


    config = PresolveOptions{T}()

    # Round the bounds of integer variables are integers.
    round_integer_bounds!(ps)

    # Strengthen the bounds of integer variables by domain propagation.
    bound_strengthening!(ps)

    # Check bound consistency on all rows/columns
    st = bounds_consistency_checks!(ps)
    ps.status == PRIMAL_INFEASIBLE && return ps.status

    # I. Remove all fixed variables, empty rows and columns
    @_return_if_inferred apply!(ps, RemoveFixedVariables(), config)
    @_return_if_inferred apply!(ps, RemoveEmptyRows(), config)
    @_return_if_inferred remove_empty_columns!(ps)

    # Identify row singletons
    ps.row_singletons = [i for (i, nz) in enumerate(ps.nzrow) if ps.rowflag[i] && nz == 1]

    # II. Passes
    ps.updated = true
    npasses = 0  # TODO: maximum number of passes
    while ps.updated && ps.status == NOT_INFERRED
        npasses += 1
        ps.updated = false
        @debug "Presolve pass $npasses" ps.nrow ps.ncol
        round_integer_bounds!(ps)

        bound_strengthening!(ps)

        @_return_if_inferred bounds_consistency_checks!(ps)
        @_return_if_inferred remove_empty_columns!(ps)


        # Remove all fixed variables
        # TODO: remove empty variables as well
        @_return_if_inferred remove_row_singletons!(ps)
        @_return_if_inferred apply!(ps, RemoveFixedVariables(), config)

        # Remove forcing & dominated constraints
        @_return_if_inferred remove_row_singletons!(ps)
        @_return_if_inferred remove_forcing_rows!(ps)

        # Remove free and implied free column singletons
        @_return_if_inferred remove_row_singletons!(ps)
        @_return_if_inferred remove_free_column_singletons!(ps)

        # TODO: remove column singleton with doubleton equation

        # Dual reductions
        @_return_if_inferred remove_row_singletons!(ps)
        @_return_if_inferred remove_dominated_columns!(ps)
    end

    @_return_if_inferred remove_empty_columns!(ps)

    @debug("Presolved problem info",
        ps.pb0.ncon, ps.nrow,
        ps.pb0.nvar, ps.ncol,
        sum(ps.nzcol[ps.colflag]), sum(ps.nzrow[ps.rowflag])
    )

    # TODO: check problem dimensions and declare optimality if problem is empty
    if ps.nrow == 0 && ps.ncol == 0
        # Problem is empty: declare optimality now
        ps.status = OPTIMAL

        # Resize solution
        resize!(ps.solution, 0, 0)
        ps.solution.primal_status = FEASIBLE_POINT
        ps.solution.dual_status = FEASIBLE_POINT
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = false
        ps.solution.z_primal = ps.obj0
        ps.solution.z_dual = ps.obj0
    end

    # Old <-> new index mapping
    compute_index_mapping!(ps)

    # TODO: extract reduced problem (?)

    # Done.
    return ps.status
end

function compute_index_mapping!(ps::PresolveData)
    ps.new_con_idx = Vector{Int}(undef, ps.pb0.ncon)
    ps.new_var_idx = Vector{Int}(undef, ps.pb0.nvar)
    ps.old_con_idx = Vector{Int}(undef, ps.nrow)
    ps.old_var_idx = Vector{Int}(undef, ps.ncol)

    inew = 0
    for iold in 1:ps.pb0.ncon
        if ps.rowflag[iold]
            inew += 1
            ps.new_con_idx[iold] = inew
            ps.old_con_idx[inew] = iold
        else
            ps.new_con_idx[iold] = 0
        end
    end
    jnew = 0
    for jold in 1:ps.pb0.nvar
        if ps.colflag[jold]
            jnew += 1
            ps.new_var_idx[jold] = jnew
            ps.old_var_idx[jnew] = jold
        else
            ps.new_var_idx[jold] = 0
        end
    end

    return nothing
end

"""
    bounds_consistency_checks(ps)

Check that all primal & dual bounds are consistent.

TODO: If not, declare primal/dual infeasibility and extract ray.
"""
function bounds_consistency_checks!(ps::PresolveData{T}) where {T}
    # Check primal bounds
    for (i, (l, u)) in enumerate(zip(ps.lrow, ps.urow))
        if ps.rowflag[i] && l > u
            # Problem is primal infeasible
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

            # Farkas ray: y⁺_i = y⁻_i = 1 (any > 0 value works)
            ps.solution.primal_status = NO_SOLUTION
            ps.solution.dual_status = INFEASIBILITY_CERTIFICATE
            ps.solution.is_primal_ray = false
            ps.solution.is_dual_ray = true
            ps.solution.z_primal = ps.solution.z_dual = T(Inf)
            i_ = ps.new_con_idx[i]
            ps.solution.y_lower[i_] = one(T)
            ps.solution.y_upper[i_] = one(T)

            return
        end
    end
    for (j, (l, u)) in enumerate(zip(ps.lcol, ps.ucol))
        if ps.colflag[j] && l > u
            # Primal is primal infeasible
            @debug "Column $j is primal infeasible"
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

            # Farkas ray: y⁺_i = y⁻_i = 1 (any > 0 value works)
            ps.solution.primal_status = NO_SOLUTION
            ps.solution.dual_status = INFEASIBILITY_CERTIFICATE
            ps.solution.is_primal_ray = false
            ps.solution.is_dual_ray = true
            ps.solution.z_primal = ps.solution.z_dual = T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.s_lower[j_] = one(T)
            ps.solution.s_upper[j_] = one(T)

            return
        end
    end

    # TODO: Check dual bounds

    return nothing
end

"""
    remove_empty_columns!(ps::PresolveData)

Remove all empty columns.

Called once at the beginning of the presolve procedure.
If an empty column is created later, it is removed on the spot.
"""
function remove_empty_columns!(ps::PresolveData{T}) where {T}
    for j in 1:ps.pb0.nvar
        remove_empty_column!(ps, j)
        ps.status == NOT_INFERRED || break
    end
    return nothing
end


"""
    round_integer_bounds!(ps::PresolveData)

Ensure all columns with integer variables having integer bounds.

Called once at the very beginning of the presolve procedure.
"""

function round_integer_bounds!(ps::PresolveData{T}) where {T}
    # The problem is LP.
    ps.pb0.is_continuous && return nothing

    for j in 1:ps.pb0.nvar
        round_integer_bounds!(ps, j)
    end
    return nothing
end

"""
    bound_strengthening!(ps::PresolveData)

Strengthen the bounds on integer variables with domain propagation.
"""

function bound_strengthening!(ps::PresolveData{T}) where {T}
    # The problem is LP.
    ps.pb0.is_continuous && return nothing

    for j in 1:ps.pb0.nvar
        bound_strengthening!(ps, j)
    end
    return nothing
end

function remove_row_singletons!(ps::PresolveData{T}) where {T}
    nsingletons = 0
    for i in ps.row_singletons
        remove_row_singleton!(ps, i)
    end
    ps.row_singletons = Int[]
    return nothing
end

"""
    remove_forcing_rows!

Remove forcing and dominated row
"""
function remove_forcing_rows!(ps::PresolveData)
    for (i, flag) in enumerate(ps.rowflag)
        flag && remove_forcing_row!(ps, i)
    end
    return nothing
end

"""
    remove_free_column_singletons!(ps)

"""
function remove_free_column_singletons!(ps::PresolveData)
    for (j, flag) in enumerate(ps.colflag)
        remove_free_column_singleton!(ps, j)
    end
    return nothing
end

function remove_dominated_columns!(ps::PresolveData{T}) where {T}
    # Strengthen dual bounds with column singletons
    for (j, (l, u)) in enumerate(zip(ps.lcol, ps.ucol))
        (ps.colflag[j] && ps.nzcol[j] == 1) || continue

        col = ps.pb0.acols[j]
        # Find non-zero index
        nz = 0
        i, aij = 0, zero(T)
        for (i_, a_) in zip(col.nzind, col.nzval)
            if ps.rowflag[i_] && !iszero(a_)
                nz += 1; nz <= 1 || break
                i = i_
                aij = a_
            end
        end

        nz == 1 || (@error "Expected singleton but column $j has $nz non-zeros"; continue)
        iszero(aij) && continue  # empty column

        # Strengthen dual bounds
        #=

        =#
        cj = ps.obj[j]
        y_ = cj / aij
        if !isfinite(l) && !isfinite(u)
            # Free variable. Should not happen but handle anyway
            # TODO
        elseif isfinite(l) && !isfinite(u)
            # Lower-bounded variable: `aij * yi ≤ cj`
            if aij > zero(T)
                # yi ≤ cj / aij
                @debug "Col $j forces y$i <= $y_"
                ps.uy[i] = min(ps.uy[i],  y_)
            else
                # yi ≥ cj / aij
                @debug "Col $j forces y$i >= $y_"
                ps.ly[i] = max(ps.ly[i],  y_)
            end

        elseif !isfinite(l) && isfinite(u)
            # Upper-bounded variable: `aij * yi ≥ cj`
            if aij > zero(T)
                # yi ≥ cj / aij
                @debug "Col $j forces y$i >= $y_"
                ps.ly[i] = max(ps.ly[i],  y_)
            else
                # yi ≤ cj / aij
                @debug "Col $j forces y$i <= $y_"
                ps.uy[i] = min(ps.uy[i],  y_)
            end
        end

        # TODO: dual feasibility check (?)
    end

    for (j, flag) in enumerate(ps.colflag)
        remove_dominated_column!(ps, j)
        ps.status == NOT_INFERRED || break
    end
    return nothing
end
