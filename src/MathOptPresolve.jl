module MathOptPresolve

using LinearAlgebra, SparseArrays

using MathOptInterface
const MOI = MathOptInterface

include("status.jl")
include("options.jl")
include("problem_data.jl")
include("solution.jl")
include("presolve.jl")
include("util.jl")

include("moi.jl")

end
