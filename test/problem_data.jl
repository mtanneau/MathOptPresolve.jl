function run_tests_pbdata(T::Type)

    @testset "Creation" begin
        pb = MathOptPresolve.ProblemData{T}("test-empty")
        @test pb.name == "test-empty"
        check_problem_size(pb, 0, 0)
        @test pb.objsense
        @test iszero(pb.obj0)

        # The model is:
        #=
            min     x1 + 2 x2
            s.t.    -∞ ⩽  -x1 +   x2 ⩽ 1
                    -1 ⩽ 2 x1 - 2 x2 ⩽ 0
                    0 ⩽ x1 ⩽ ∞
                    1 ⩽ x2 ⩽ ∞
                    x2 in \{0,1\} =#

        MathOptPresolve.load_problem!(
            pb,
            "test-filled",
            true,
            T[1.0, 2.0],
            zero(T),
            sparse([-1.0  1.0
                     2.0 -2.0]),
            T[-Inf, -1.0],
            T[1.0, 0.0],
            T[0.0, 1.0],
            T[Inf, Inf],
            MathOptPresolve.VariableType[MathOptPresolve.CONTINUOUS, MathOptPresolve.BINARY],
        )

        @test pb.name == "test-filled"

        @test pb.objsense
        @test iszero(pb.obj0)

        col1, col2 = pb.acols[1], pb.acols[2]
        @test pb.obj == [one(T), 2 * one(T)]
        @test pb.lvar == [zero(T), one(T)]
        @test pb.uvar == [T(Inf), T(Inf)]

        # Check dimensions
        check_problem_size(pb, 2, 2)

        # Check coefficients
        row1, row2 = pb.arows[1], pb.arows[2]
        @test row1.nzind == [1, 2]
        @test row1.nzval == T.([-1, 1])
        @test row2.nzind == [1, 2]
        @test row2.nzval == T.([2, -2])
        @test col1.nzind == [1, 2]
        @test col1.nzval == T.([-1, 2])
        @test col2.nzind == [1, 2]
        @test col2.nzval == T.([1, -2])

        # Check row bounds
        @test pb.lcon == [T(-Inf), -one(T)]
        @test pb.ucon == [one(T), zero(T)]

        @test pb.var_types == [MathOptPresolve.CONTINUOUS, MathOptPresolve.BINARY]
        @test !pb.is_continuous

        empty!(pb)
        @test pb.name == ""
        @test iszero(pb.obj0)
        check_problem_size(pb, 0, 0)
    end

    @testset "presolve!" begin
        pb = MathOptPresolve.ProblemData{T}()

        # The model is:
        #=
            min     x1 + 2 x2
            s.t.    1 ⩽  x1 ⩽ 1
                    0 ⩽ x1 ⩽ ∞
                    1 ⩽ x2 ⩽ ∞ =#

        MOP.load_problem!(
            pb,
            "optimal",
            true,
            T[1.0, 2.0],
            zero(T),
            sparse([1.0  0.0]),
            T[1.0],
            T[1.0],
            T[0.0, 1.0],
            T[Inf, Inf],
        )
        ps = MOP.PresolveData(pb)
        status = MOP.presolve!(ps)

        @test status == MOP.OPTIMAL
        @test ps.solution.m == 0
        @test ps.solution.n == 0
        @test ps.solution.primal_status == MOP.FEASIBLE_POINT
        @test ps.solution.dual_status == MOP.FEASIBLE_POINT
        @test ps.solution.z_primal == 3.0
        @test ps.solution.z_dual == 3.0
    end

    return nothing
end

function check_problem_size(pb::MathOptPresolve.ProblemData, ncon::Int, nvar::Int)
    @test pb.ncon == ncon
    @test pb.nvar == nvar

    @test length(pb.obj) == nvar

    @test length(pb.arows) == ncon
    @test length(pb.acols) == nvar

    @test length(pb.lcon) == ncon
    @test length(pb.ucon) == ncon
    @test length(pb.lvar) == nvar
    @test length(pb.uvar) == nvar

    return nothing
end

@testset "ProblemData" begin
    for T in COEFF_TYPES
        @testset "$T" begin run_tests_pbdata(T) end
    end
end
