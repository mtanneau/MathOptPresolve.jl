function _postsolve_fn(ps::PresolveData)
    return (orig_sol, transformed_sol) -> postsolve!(ps, orig_sol, transformed_sol)
end

mutable struct IntermediaryRepresentation{T<:Real}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
    lcon::Vector{T}
    ucon::Vector{T}
    lvar::Vector{T}
    uvar::Vector{T}
    var_types::Vector{VariableType}

    function IntermediaryRepresentation{T}(n::Int) where {T<:Real}
        return new{T}(
            Int[],
            Int[],
            T[],
            T[],
            T[],
            fill(T(-Inf), n),
            fill(T(Inf), n),
            fill(CONTINUOUS, n),
        )
    end
end

num_rows(ir::IntermediaryRepresentation) = length(ir.lcon)

function add_row!(ir::IntermediaryRepresentation)::Int
    if isempty(ir.I)
        @assert isempty(ir.J)
        @assert isempty(ir.V)
        @assert isempty(ir.lcon)
        @assert isempty(ir.ucon)
        return 1
    else
        row_num = I[end]
        @assert row_num == length(lcon) == length(ucon)
        return row_num + 1
    end
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    s::MOI.EqualTo{T},
) where {T<:Real}
    i = f.variable.value
    ir.lvar[i] = s.value
    ir.uvar[i] = s.value
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    s::MOI.LessThan{T},
) where {T<:Real}
    i = f.variable.value
    ir.uvar[i] = s.upper
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    s::MOI.GreaterThan{T},
) where {T<:Real}
    i = f.variable.value
    ir.lvar[i] = s.lower
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    s::MOI.Interval{T},
) where {T<:Real}
    i = f.variable.value
    ir.lvar[i] = s.lower
    ir.uvar[i] = s.upper
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    ::MOI.ZeroOne,
)
    i = f.variable.value
    ir.var_types[i] = BINARY
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.SingleVariable,
    ::MOI.Integer,
)
    i = f.variable.value
    ir.var_types[i] = GENERAL_INTEGER
    return nothing
end

_get_bounds(s::MOI.EqualTo) = (s.value, s.value)
_get_bounds(s::MOI.LessThan{T}) where {T} = (T(-Inf), s.upper)
_get_bounds(s::MOI.GreaterThan{T}) where {T} = (s.lower, T(Inf))
_get_bounds(s::MOI.Interval) = (s.lower, s.upper)
_get_bounds(s) = error("Unexpected set $s.")

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.ScalarAffineFunction{T},
    s::MOI.AbstractScalarSet,
) where {T<:Real}
    row_id = add_row!(ir)
    for term in f.terms
        vi = term.variable_index
        push!(ir.I, row_id)
        push!(ir.J, term.variable_index.value)
        push!(ir.V, term.coefficient)
    end
    lb, ub = _get_bounds(s)
    push!(ir.lcon, lb - f.constant)
    push!(ir.ucon, ub - f.constant)
    return nothing
end

process_constraint!(ir::IntermediaryRepresentation, f, s) =
    error("Unsupported constraint of type $(typeof(f))-in-$(typeof(s)).")

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.AbstractVectorFunction,
    s::MOI.AbstractVectorSet,
)
    return process_constraint!(ir, MOIU.scalarize(f), MOIU.scalarize(s))
end

function presolve!(dest::MOI.ModelLike, src::MOI.ModelLike, T::Type{<:Real})
    @assert MOI.is_empty(dest)

    objsense = MOI.get(src, MOI.ObjectiveSense())
    obj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())

    n = MOI.get(src, MOI.NumberOfVariables())
    ir = IntermediaryRepresentation{T}(n)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(src, MOI.ConstraintFunction(), ci)
            s = MOI.get(src, MOI.ConstraintSet(), ci)
            process_constraint!(ir, f, s)
        end
    end
    obj_vec = zeros(T, n)
    for term in obj.terms
        obj_vec[term.variable_index.value] += term.coefficient
    end
    pd = ProblemData{T}()
    load_problem!(
        pd,
        "model-from-MOI",
        objsense == MOI.MIN_SENSE,
        obj_vec,
        MOI.constant(obj),
        SparseArrays.dropzeros!(SparseArrays.sparse(ir.I, ir.J, ir.V, num_rows(ir), n)),
        ir.lcon,
        ir.ucon,
        ir.lvar,
        ir.uvar,
        ir.var_types,
    )
    ps = PresolveData(pd)
    presolve!(ps)
    MOI.empty!(dest)
    if ps.status == OPTIMAL
        return MOI.OPTIMAL, x -> _optimal_postsolve_fn(ps, x)
    elseif ps.status == PRIMAL_INFEASIBLE
        return MOI.PRIMAL_INFEASIBLE, x -> _infeasible_postsolve_fn(ps, x)
    elseif ps.status == DUAL_INFEASIBLE
        return MOI.DUAL_INFEASIBLE, x -> _unbounded_postsolve_fn(ps, x)
    else
        @assert ps.status == NOT_INFERRED
        @assert moi_status == MOI.OPTIMIZE_NOT_CALLED
        extract_reduced_problem!(ps)

        pd = ps.pb_red
        @assert pd.nvar > 0
        MOI.empty!(dest)

        vis = MOI.add_variables(dest, pd.nvar)
        x = [MOI.SingleVariable(vis[j]) for j = 1:pd.nvar]
        for j = 1:pd.nvar
            MOI.add_constraint(dest, x[j], MOI.Interval{T}(pd.lvar[j], pd.lcon[j]))
            if pd.var_types[j] == BINARY
                MOI.add_constraint(dest, x[j], MOI.ZeroOne())
            elseif pd.var_types[j] == GENERAL_INTEGER
                MOI.add_constraint(dest, x[j], MOI.Integer())
            end
        end
        MOI.set(dest, MOI.ObjectiveSense(), pd.objsense ? MOI.MIN_SENSE : MOI.MAX_SENSE)
        MOI.set(
            dest,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
            sum(pd.obj[i] * x[i] for i = 1:pd.nvar) + pd.obj0,
        )
        for i = 1:pd.ncon
            row = pd.arows[i]
            lb, ub = pd.lcon[i], pd.ucon[i]
            set = (
                if lb == T(-Inf)
                    MOI.LessThan{T}(ub)
                elseif ub == T(Inf)
                    MOI.GreaterThan{T}(lb)
                elseif lb ≈ ub
                    MOI.EqualTo{T}(lb)
                else
                    MOI.Interval{T}(lb, ub)
                end
            )
            MOI.add_constraint(
                dest,
                sum(row.nzval[k] * x[row.nzind[k]] for k = 1:length(row.nzind)),
                scalar_set,
            )
        end
        return MOI.OPTIMIZE_NOT_CALLED, x -> _not_inferred_postsolve_fn(ps, x)
    end
end

function _optimal_postsolve_fn(ps::PresolveData{T}, x::Vector{T}) where {T}
    if !isempty(x)
        throw(ArgumentError("Presolve solved model to optimality; postsolve expects an empty input argument."))
    end
    @assert ps.solution.primal_status == FEASIBLE_POINT
    @assert !ps.solution.is_primal_ray
    @assert !ps.solution.is_dual_ray
    return ps.solution.x
end

function _infeasible_postsolve_fn(ps::PresolveData{T}, x::Vector{T}) where {T}
    if !isempty(x)
        throw(ArgumentError("Presolve proved that model is infeasible; postsolve expects an empty input argument."))
    end
    @assert ps.solution.primal_status == NO_SOLUTION
    @assert !ps.solution.is_primal_ray
    if ps.solution.is_dual_ray
        # TODO: Return dual ray here
        return nothing
    else
        return nothing
    end
end

function _unbounded_postsolve_fn(ps::PresolveData{T}, x::Vector{T}) where {T}
    if !isempty(x)
        throw(ArgumentError("Presolve proved that model is unbounded; postsolve expects an empty input argument."))
    end
    @assert ps.solution.dual_status == NO_SOLUTION
    @assert !ps.solution.is_dual_ray
    if ps.solution.is_primal_ray
        return ps.solution.x
    else
        return nothing
    end
end

function _not_inferred_postsolve_fn(ps::PresolveData{T}, x::Vector{T}) where {T}
    orig_sol = Solution(ps.pb0.m, ps.pb0.n)
    trans_sol = Solution(ps.nrow, ps.ncol)
    trans_sol.primal_status = FEASIBLE_POINT
    trans_sol.x = x
    postsolve!(orig_sol, trans_sol, ps)
    @assert orig_sol.primal_status == FEASIBLE_POINT
    @assert !ps.solution.is_primal_ray
    @assert !ps.solution.is_dual_ray
    return orig_sol.x
end
