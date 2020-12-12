"""
    ModelStatus

- `NOT_INFERRED`: Presolve has not inferred the status of the problem.
- `OPTIMAL`: The current solution(s) are optimal.
- `PRIMAL_INFEASIBLE`: The model is primal infeasible.
- `DUAL_INFEASIBLE`: The model is dual infeasible.
"""
@enum(ModelStatus,
    NOT_INFERRED,
    OPTIMAL,
    PRIMAL_INFEASIBLE,
    DUAL_INFEASIBLE,
)

"""
    SolutionStatus

Solution Status code
- `NO_SOLUTION`: No current point.
- `FEASIBLE_POINT`: The current point is feasible.
- `INFEASIBLE_POINT`: The current point is infeasible.
- `INFEASIBILITY_CERTIFICATE`: The current point is a proof of infeasibility
of the other problem. For example, a dual solution with INFEASIBILITY_CERTIFICATE
status is a proof of infeasibility of the primal.
"""
@enum(SolutionStatus,
    NO_SOLUTION,
    FEASIBLE_POINT,
    INFEASIBLE_POINT,
    INFEASIBILITY_CERTIFICATE,
)
