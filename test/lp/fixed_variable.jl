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
Fix an integer variable to a fractional value.
"""
function test_fix_integer_fractional(T::Type)

    pb = MathOptPresolve.ProblemData{T}()

    m, n = 0, 1
    arows = Int[]
    acols = Int[]
    avals = T[]
    A = sparse(arows, acols, avals, m, n)

    # Binary variable
    MathOptPresolve.load_problem!(pb, "Test",
        true, zeros(T, n), zero(T),
        A,
        zeros(T, m), ones(T, m),
        [T(1 // 2)], [T(1 // 2)]
    )

    pb.var_types[1] = MOP.CONTINUOUS
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.NOT_INFERRED
    @test !ps.colflag[1]

    pb.var_types[1] = MOP.BINARY
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.PRIMAL_INFEASIBLE

    pb.var_types[1] = MOP.GENERAL_INTEGER
    ps = MathOptPresolve.PresolveData(pb)
    MathOptPresolve.remove_fixed_variable!(ps, 1)
    @test ps.status == MOP.PRIMAL_INFEASIBLE

    return nothing
end

@testset "Fixed variable" begin
    for T in COEFF_TYPES
        @testset "$T" begin
            test_fixed_variable_with_zeros(T)
            test_fix_integer_fractional(T)
        end
    end
end
