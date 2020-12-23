function round_integer_bounds_tests(T::Type)
    # Build the following model
    #=
        min     x + y + z + w
        s.t.   -10 ⩽ 2 x + 3 y + 4 z + w ⩽ 10
                -1.5 ⩽ x ⩽ 1.6
                -2.5 ⩽ y ⩽ 2.1
                -1.3 ⩽ z ⩽ 3
                MOP.eps(T)/2 ⩽ w ⩽ 1 - MOP.eps(T)/2
                x is integer, z, w are binary  =#
    C = T[1, 1, 1, 1]
    lc = T[-3 // 2, -5 // 2, -13 // 10, MOP.eps(T)/2]
    uc = T[8 // 5, 21 // 10, 3, 1 - MOP.eps(T)/2]
    lr = T[-10]
    ur = T[10]
    varTypes = [MOP.GENERAL_INTEGER, MOP.CONTINUOUS,
                MOP.BINARY, MOP.BINARY]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse([2 3 4 1]), lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)

    MOP.round_integer_bounds!(ps)

    @test ps.lcol == T.([-1, -5 // 2, 0, 0])
    @test ps.ucol == T.([1, 21 // 10, 1, 1])
    @test ps.var_types == [MOP.GENERAL_INTEGER, MOP.CONTINUOUS,
                MOP.BINARY, MOP.BINARY]
    @test ps.pb0.var_types == [MOP.GENERAL_INTEGER, MOP.CONTINUOUS,
                MOP.BINARY, MOP.BINARY]
    @test ps.status == MOP.NOT_INFERRED
    return nothing
end

function round_integer_bounds_tests_2(T::Type)
    # Build the following model
    #=
        min     x + y
        s.t.    MOP.eps(T)/2 ⩽ x ⩽ MOP.eps(T)/2
                0.9 ⩽ y ⩽ 0.1
                x, y are integer  =#
    C = T.([1, 1])
    lc = T.([MOP.eps(T)/2, 0.9])
    uc = T.([MOP.eps(T)/2, 0.1])
    lr = Vector{T}([])
    ur = Vector{T}([])
    varTypes = [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        spzeros(T, 0, 2), lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)

    MOP.round_integer_bounds!(ps)

    @test ps.lcol == T.([0, 1])
    @test ps.ucol == T.([0, 0])
    @test ps.pb0.var_types == [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER]
    @test ps.var_types == [MOP.BINARY, MOP.BINARY]
    @test ps.status == MOP.NOT_INFERRED
    return nothing
end

@testset "Round Integer Bounds" begin
    for T in COEFF_TYPES
        @testset "$T" begin round_integer_bounds_tests(T) end
        @testset "$T" begin round_integer_bounds_tests_2(T) end
    end
end
