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
        row_num = ir.I[end]
        @assert row_num == length(ir.lcon) == length(ir.ucon)
        return row_num + 1
    end
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    s::MOI.EqualTo{T},
) where {T<:Real}
    i = f.value
    ir.lvar[i] = s.value
    ir.uvar[i] = s.value
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    s::MOI.LessThan{T},
) where {T<:Real}
    i = f.value
    ir.uvar[i] = s.upper
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    s::MOI.GreaterThan{T},
) where {T<:Real}
    i = f.value
    ir.lvar[i] = s.lower
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    s::MOI.Interval{T},
) where {T<:Real}
    i = f.value
    ir.lvar[i] = s.lower
    ir.uvar[i] = s.upper
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    ::MOI.ZeroOne,
)
    i = f.value
    ir.var_types[i] = BINARY
    return nothing
end

function process_constraint!(
    ir::IntermediaryRepresentation,
    f::MOI.VariableIndex,
    ::MOI.Integer,
)
    i = f.value
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
        push!(ir.I, row_id)
        push!(ir.J, term.variable.value)
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

struct PresolveResult{T}
    pd::PresolveData{T}
end

function get_status(pr::PresolveResult)
    if pr.pd.status == OPTIMAL
        return MOI.OPTIMAL
    elseif pr.pd.status == PRIMAL_INFEASIBLE
        return MOI.INFEASIBLE
    elseif pr.pd.status == DUAL_INFEASIBLE
        return MOI.DUAL_INFEASIBLE
    else
        @assert pr.pd.status == NOT_INFERRED
        return MOI.OPTIMIZE_NOT_CALLED
    end
end

"""
Query the optimal solution produced by presolve. The solution returned is for
the original, not transformed, problem.

Note: Throws an error if presolve did not terminate with status `MOI.OPTIMAL`.
"""
function get_optimal_solution(pr::PresolveResult{T}) where {T}
    if get_status(pr) != MOI.OPTIMAL
        error(
            "Presolve did not solve problem to optimality, so you cannot query an optimal solution.",
        )
    end
    @assert pr.pd.solution !== nothing
    @assert pr.pd.solution.primal_status == FEASIBLE_POINT
    @assert pr.pd.nrow == pr.pd.ncol == 0
    orig_sol = _post_crush(pr.pd, pr.pd.solution)
    @assert orig_sol.primal_status == FEASIBLE_POINT
    return orig_sol.x
end

"""
Query the unbounded ray produced by presolve. The ray returned is for the
original, not transformed, problem.

Note: Throws an error if presolve did not terminate with status
`MOI.DUAL_INFEASIBLE`.
"""
function get_unbounded_ray(pr::PresolveResult{T}) where {T}
    if get_status(pr) != MOI.DUAL_INFEASIBLE
        error("Presolve did not prove unboundedness, so you cannot query an unbounded ray.")
    end
    @assert pr.pd.solution !== nothing
    @assert pr.pd.solution.is_primal_ray
    @assert pr.pd.solution.primal_status == INFEASIBILITY_CERTIFICATE
    orig_sol = _post_crush(pr.pd, pr.pd.solution)
    @assert orig_sol.primal_status == INFEASIBILITY_CERTIFICATE
    return orig_sol.x
end

function get_infeasibility_certificate(pr::PresolveResult{T}) where {T}
    if get_status(pr) != MOI.INFEASIBLE
        error(
            "Presolve did not prove infeasibility, so you cannot query an infeasibility certificate.",
        )
    end
    error("Not yet implemented. Come back soon!")
end

"""
    post_crush(pr::PresolveResult{T}, x::Vector{T}; is_ray::Bool=false)

Map a point from the presolved model to the original model. If
`is_ray = false`, the point should be a feasible point for the presolved model;
otherwise, it should be an unbounded ray.

Note: Throws an error if the length of `x` is not equal to the number of
variables in the presolved model.
"""
function post_crush(pr::PresolveResult{T}, x::Vector{T}; is_ray::Bool = false) where {T}
    if length(x) != pr.pd.ncol
        throw(
            ArgumentError(
                "Transformed solution is of length $(length(x)); expected one of length $(pr.pd.ncol).",
            ),
        )
    end
    trans_sol = Solution{T}(pr.pd.nrow, pr.pd.ncol)
    trans_sol.primal_status = is_ray ? INFEASIBILITY_CERTIFICATE : FEASIBLE_POINT
    trans_sol.x = x
    trans_sol.is_primal_ray = is_ray
    return _post_crush(pr.pd, trans_sol).x
end

function _post_crush(pd::PresolveData{T}, trans_sol::Solution{T}) where {T}
    orig_sol = Solution{T}(pd.pb0.ncon, pd.pb0.nvar)
    postsolve!(orig_sol, trans_sol, pd)
    return orig_sol
end

"""
    presolve!(dest::MOI.ModelLike, src::MOI.ModelLike, T::Type{<:Real})::PresolveResult

Apply presolve to `src` model, and populate `dest` with the reduced problem.
The type `T` specifies the data type used to represent the objective,
constraints, etc. added to `dest`.

Returns a `PresolveResult` result object. You may query the status of the model
after the presolve routine (`get_status`), and depending on this status, also
query the relevant solution information (with `get_optimal_solution`,
`get_unbounded_ray`, and `get_infeasibility_certificate`).

Note: The function will throw an `ArgumentError` if `dest` is not empty.
"""
function presolve!(dest::MOI.ModelLike, src::MOI.ModelLike, T::Type{<:Real})
    if !MOI.is_empty(dest)
        throw(ArgumentError("Destination model is not empty."))
    end

    objsense = MOI.get(src, MOI.ObjectiveSense())
    obj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())

    n = MOI.get(src, MOI.NumberOfVariables())
    ir = IntermediaryRepresentation{T}(n)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(src, MOI.ConstraintFunction(), ci)
            s = MOI.get(src, MOI.ConstraintSet(), ci)
            process_constraint!(ir, f, s)
        end
    end
    obj_vec = zeros(T, n)
    for term in obj.terms
        obj_vec[term.variable.value] += term.coefficient
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
    if ps.status != OPTIMAL
        extract_reduced_problem!(ps)
        pd = ps.pb_red
        vis = MOI.add_variables(dest, pd.nvar)
        x = [vis[j] for j = 1:pd.nvar]
        for j = 1:pd.nvar
            MOI.add_constraint(dest, x[j], MOI.Interval{T}(pd.lvar[j], pd.uvar[j]))
            if pd.var_types[j] == BINARY
                MOI.add_constraint(dest, x[j], MOI.ZeroOne())
            elseif pd.var_types[j] == GENERAL_INTEGER
                MOI.add_constraint(dest, x[j], MOI.Integer())
            end
        end
        MOI.set(dest, MOI.ObjectiveSense(), pd.objsense ? MOI.MIN_SENSE : MOI.MAX_SENSE)
        MOI.set(
            dest,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
            sum(pd.obj[i] * x[i] for i = 1:pd.nvar) + pd.obj0,
        )
        for i = 1:pd.ncon
            row = pd.arows[i]
            lb, ub = pd.lcon[i], pd.ucon[i]
            scalar_set = (
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
    end
    return PresolveResult{T}(ps)
end
