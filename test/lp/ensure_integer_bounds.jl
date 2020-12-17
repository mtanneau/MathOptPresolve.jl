function ensure_integer_bounds_tests(T::Type)

    # Build the following model
    #=
        min     x + y + z + w
        s.t.   -10 ⩽ 2 x + 3 y + 4 z + w ⩽ 10
                -1.5 ⩽ x ⩽ 1.6
                -2.5 ⩽ y ⩽ 2.1
                -1.3 ⩽ z ⩽ 3
                0 ⩽ w ⩽ 0
                x is integer, z, w are binary  =#
    C = T.([1, 1, 1, 1])
    lc = T.([-1.5, -2.5, -1.3, 0])
    uc = T.([1.6, 2.1, 3, 0])
    lr = T.([-10])
    ur = T.([10])
    varTypes = [MathOptPresolve.GENERAL_INTEGER, MathOptPresolve.CONTINUOUS,
                MathOptPresolve.BINARY]

    pb = MathOptPresolve.ProblemData{T}()

    MathOptPresolve.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(T.([2,3,4,1])), lc, uc, lr, ur,
        var_types=varTypes
    )

    ps = MathOptPresolve.PresolveData(pb)

    MathOptPresolve.ensure_integer_bounds(ps)

    @test ps.lcol == T([-1, -2.5, 0, 0])
    @test ps.ucol == T([1, 2.1, 1, 0])
    return
end

@testset "Ensure Integer Bounds" begin
    for T in COEFF_TYPES
        @testset "$T" begin emtpy_column_tests(T) end
    end
end
