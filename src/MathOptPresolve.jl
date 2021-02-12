module MathOptPresolve

using LinearAlgebra, SparseArrays

using MathOptInterface
const MOI = MathOptInterface

export PresolveResult,
    get_status,
    get_optimal_solution,
    get_unbounded_ray,
    get_infeasibility_certificate,
    post_crush

include("status.jl")
include("options.jl")
include("problem_data.jl")
include("solution.jl")
include("presolve.jl")
include("util.jl")

include("moi.jl")

end
