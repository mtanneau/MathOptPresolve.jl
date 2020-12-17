"""
Remove fixed variables with explicit zeros.
"""
function test_fixed_variable_with_zeros(T::Type)

    pb = MathOptPresolve.ProblemData{T}()

    m, n = 3, 2
    arows = [1, 1, 2, 2, 3, 3]
    acols = [1, 2, 1, 2, 1, 2]
    avals = T.([
        11, 0,
        0, 21,
        31, 0
    ])
    A = sparse(arows, acols, avals, m, n)

    MathOptPresolve.load_problem!(pb, "Test",
        true, T.(collect(1:n)), zero(T),
        A,
        zeros(T, m), ones(T, m),
        ones(T, n), ones(T, n),
    )

    ps = MathOptPresolve.PresolveData(pb)

    @test ps.nzrow == [1, 1, 1]
    @test ps.nzcol == [2, 1]

    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.colflag == [false, true]
    @test ps.obj0 == T(1)
    @test ps.nzrow == [0, 1, 0]

    MathOptPresolve.remove_fixed_variable!(ps, 2)
    @test ps.colflag == [false, false]
    @test ps.obj0 == T(3)
    @test ps.nzrow == [0, 0, 0]


    return nothing
end

"""
    test_fix_integer_infeasibility(T)

Test integer infeasibility when fixing integer variables.
"""
function test_fix_integer_infeasibility(T::Type)

    pb = MathOptPresolve.ProblemData{T}()

    m, n = 0, 1
    arows = Int[]
    acols = Int[]
    avals = T[]
    A = sparse(arows, acols, avals, m, n)

    MathOptPresolve.load_problem!(pb, "Test",
        true, zeros(T, n), zero(T),
        A,
        zeros(T, m), ones(T, m),
        zeros(T, n), zeros(T, n)
    )

    # Single variable with bounds ¹/₂ ≤ x ≤ ¹/₂
    pb.lvar[1] = T(1 // 2)
    pb.uvar[1] = T(1 // 2)

    # x is continuous --> OK
    pb.var_types[1] = MOP.CONTINUOUS
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.NOT_INFERRED
    @test !ps.colflag[1]

    # x in binary --> infeasible
    pb.var_types[1] = MOP.BINARY
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.PRIMAL_INFEASIBLE

    # x is integer --> infeasible
    pb.var_types[1] = MOP.GENERAL_INTEGER
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.PRIMAL_INFEASIBLE


    # Binary variable with integer bounds
    for x_ in [T(-1), T(0), T(1), T(2)]
        pb.lvar[1] = x_
        pb.uvar[1] = x_

        # x is integer --> OK
        pb.var_types[1] = MOP.GENERAL_INTEGER
        ps = MathOptPresolve.PresolveData(pb)
        MathOptPresolve.remove_fixed_variable!(ps, 1)
        @test ps.status == MOP.NOT_INFERRED
        @test !ps.colflag[1]

        # x in binary --> infeasible
        pb.var_types[1] = MOP.BINARY
        ps = MathOptPresolve.PresolveData(pb)
        MathOptPresolve.remove_fixed_variable!(ps, 1)
        if x_ in [T(0), T(1)]
            @test ps.status == MOP.NOT_INFERRED
        else
            @test ps.status == MOP.PRIMAL_INFEASIBLE
        end
    end

    return nothing
end

@testset "Fixed variable" begin
    for T in COEFF_TYPES
        @testset "$T" begin
            test_fixed_variable_with_zeros(T)
            test_fix_integer_infeasibility(T)
        end
    end
end
