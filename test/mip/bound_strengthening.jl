function bound_strengthening_Integer_tests(T::Type)

    config = MOP.PresolveOptions{T}()

    # Build the following model
    #=
        min     x + y + z
        s.t.    -10 + MOP.eps(T)/2 ⩽ x - 2 y + 0 z ⩽ 10 - MOP.eps(T)/2
                x, z are integer, y is binary  =#
    C = T.([1, 1, 1])
    lc = T.([-Inf, 0, -Inf])
    uc = T.([Inf, 1, Inf])
    lr = Vector{T}([-10 + MOP.eps(T)/2])
    ur = Vector{T}([10 - MOP.eps(T)/2])
    varTypes = [MOP.GENERAL_INTEGER, MOP.BINARY, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1], [1, 2], [1, -2], 1, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([-10, 0, -Inf])
    @test ps.ucol == T.([12, 1, Inf])


    # Build the following model
    #=
        min     x + y + z
        s.t.    -10 ⩽ 10 x - y + 0 z ⩽ 10
                x, y, z are integer  =#
    C = T.([1, 1, 1])
    lc = T.([-Inf, -Inf, -1])
    uc = T.([Inf, Inf, 1])
    lr = Vector{T}([-10])
    ur = Vector{T}([10])
    varTypes = [MOP.GENERAL_INTEGER, MOP.BINARY, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1], [1, 2], [10, -1], 1, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([-Inf, -Inf, -1])
    @test ps.ucol == T.([Inf, Inf, 1])

    # Build the following model
    #=
        min     x + y + z
        s.t.    -10 ⩽ 10 x - y + 0 z ⩽ 10
                x, y, z are integer  =#
    C = T.([1, 1, 1])
    lc = T.([-Inf, -Inf, -1])
    uc = T.([Inf, Inf, 1])
    lr = Vector{T}([-10])
    ur = Vector{T}([10])
    varTypes = [MOP.GENERAL_INTEGER, MOP.BINARY, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1], [1, 2], [10, -1], 1, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([-Inf, -Inf, -1])
    @test ps.ucol == T.([Inf, Inf, 1])

    # Build the following model
    #=
        min     x + y + z
        s.t.    0 ⩽ x + y - z ⩽ 1
                -1 ⩽ x ⩽ 1
                x, y, z are integer  =#
    C = T.([1, 1, 1])
    lc = T.([-1, -Inf, -Inf])
    uc = T.([1, Inf, Inf])
    lr = Vector{T}([0])
    ur = Vector{T}([1])
    varTypes = [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1, 1], [1, 2, 3], [1, 1, -1], 1, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([-1, -Inf, -Inf])
    @test ps.ucol == T.([1, Inf, Inf])
    # Build the following model
    #=
        min     x + y + z
        s.t.    MOP.eps(T)/2 ⩽ x + z ⩽ MOP.eps(T)/2
                -MOP.eps(T)/2 ⩽ y + z ⩽ -MOP.eps(T)/2
                MOP.eps(T)/2 ⩽ x + y ⩽ -MOP.eps(T)/2
                0 ⩽ x ⩽ 1
                x, y, z are integer  =#
    C = T.([1, 1, 1])
    lc = Vector{T}([0, -Inf, -Inf])
    uc = Vector{T}([1, Inf, Inf])
    lr = T.([MOP.eps(T)/2, -MOP.eps(T)/2, MOP.eps(T)/2])
    ur = T.([MOP.eps(T)/2, -MOP.eps(T)/2, -MOP.eps(T)/2])

    varTypes = [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1, 2, 2, 3, 3], [1, 3, 2, 3, 1, 2], [1, 1, 1, 1, 1, 1], 3, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([0, -1, 0])
    @test ps.ucol == T.([1, 0, 0])

    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([0, 0, 0])
    @test ps.ucol == T.([0, 0, 0])
    return nothing
end


function bound_strengthening_MIP_tests(T::Type)

    config = MOP.PresolveOptions{T}()

    # Build the following model
    #=
        min     x + y + z
        s.t.    -10.5 ⩽ x - 2 y + 0 z ⩽ 10.5
                x, z are integer, y is binary  =#
    C = T.([1, 1, 1])
    lc = T.([-Inf, 0, -Inf])
    uc = T.([Inf, 1, Inf])
    lr = Vector{T}([-21 // 2])
    ur = Vector{T}([21 // 2])
    varTypes = [MOP.CONTINUOUS, MOP.BINARY, MOP.GENERAL_INTEGER]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1], [1, 2], [1, -2], 1, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([-21 // 2, 0, -Inf])
    @test ps.ucol == T.([25 // 2, 1, Inf])

    # Build the following model
    #=
        min     x + y + z
        s.t.    MOP.eps(T)/2 ⩽ x + z ⩽ MOP.eps(T)/2
                -MOP.eps(T)/2 ⩽ y + z ⩽ -MOP.eps(T)/2
                MOP.eps(T)/2 ⩽ x + y ⩽ -MOP.eps(T)/2
                0 ⩽ x ⩽ 1
                x, y, z are integer  =#
    C = T.([1, 1, 1])
    lc = Vector{T}([0, -Inf, -Inf])
    uc = Vector{T}([1, Inf, Inf])
    lr = T.([0, 0, 0])
    ur = T.([0, 0, 0])

    varTypes = [MOP.CONTINUOUS, MOP.GENERAL_INTEGER, MOP.CONTINUOUS]
    pb = MOP.ProblemData{T}()
    A = sparse([1, 1, 2, 2, 3, 3], [1, 3, 2, 3, 1, 2], [1, 1, 1, 1, 1, 1], 3, 3)
    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        A, lr, ur, lc, uc,
        varTypes
    )

    ps = MOP.PresolveData(pb)
    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([0, -1, 0])
    @test ps.ucol == T.([1, 0, 0])

    MOP.apply!(ps, MOP.StrengthenBounds(), config)
    @test ps.lcol == T.([0, 0, 0])
    @test ps.ucol == T.([0, 0, 0])
    return nothing
end

@testset "Bound Strengthening" begin
    for T in COEFF_TYPES
        @testset "$T" begin bound_strengthening_Integer_tests(T) end
        @testset "$T" begin bound_strengthening_MIP_tests(T) end
    end
end
