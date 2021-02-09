using MathOptPresolve
using Test

using SparseArrays

const MOP = MathOptPresolve

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

const COEFF_TYPES = [Float64, BigFloat]

@testset "Problem data" begin
    include("problem_data.jl")
end

@testset "LP" begin
    @testset "empty column" begin
        include("lp/empty_column.jl")
    end
    @testset "empty row" begin
        include("lp/empty_row.jl")
    end
    @testset "fixed variable" begin
        include("lp/fixed_variable.jl")
    end
end

@testset "MIP" begin
    @testset "round integer bounds" begin
        include("mip/round_integer_bounds.jl")
    end
    @testset "coefficient strengthening" begin
        include("mip/coefficient_strengthening.jl")
    end
end

@testset "MOI" begin
    include("moi.jl")
end
