function cg_strengthening_test1(T::Type)
    # We test cg_strengthening! function on the following MIP
    #=
    min     x₁ + x₂ + x₃ + x₄ + x₅ + x₆
    s.t.    x₁ + x₂ + x₃                 ⩾ 3.7
            2x₁    + 2x₃ + 4x₄           ⩾ 8.5
            x₁ + x₂ + 3x₃ + x₄ + 7x₅     ⩾ 6.8
            -x₁ + x₂ - x₃ + x₄           ⩾ 0
            x₁                    + x₆   ⩾ 2.3

            x₁,x₂,x₃,x₄,x₅ ⩾ 0
            x₆ ⩽ 10

            x₁,x₂,x₃,x₄,x₆ are integers
            x₅ is continuous
    =#

    #= Expected outcome
    x₁ + x₂ + x₃             ⩾ 3.7 --> 4x₁ + 4x₂ + 4x₃          ⩾ 15
    2x₁    + 2x₃ + 4x₄       ⩾ 8.5 --> x₁    + x₃ + 2x₄         ⩾ 5
    We will skip the the last 3 constraints because x₅ is continous,
    x₁ is not bounded above and x₆ is not bounded below
    =#

    C = T[1, 1, 1, 1, 1, 1]
    lc = T[0, 0, 0, 0, 0, -Inf]
    uc = T[Inf, Inf, Inf, Inf, Inf, 10]
    lr = T[37//10, 85//10, 68//10, 0, 23//10]
    ur = T[Inf, Inf, Inf, Inf, Inf]
    A = T[1 1 1 0 0 0
          2 0 2 4 0 0
          1 1 3 1 7 0
          -1 1 -1 1 0 0
          1 0 0 0 0 1]

    varTypes = [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER,
                MOP.GENERAL_INTEGER, MOP.CONTINUOUS, MOP.GENERAL_INTEGER]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(A), lr, ur, lc, uc,
        varTypes
    )
    ps = MOP.PresolveData(pb)

    MOP.cg_strengthening!(ps)

    @test ps.urow == T[Inf, Inf, Inf, Inf, Inf]
    @test ps.lrow == T[15, 5, 68//10, 0, 23//10]

    @test ps.pb0.arows[1].nzval == T[4, 4, 4]
    @test ps.pb0.arows[2].nzval == T[1, 1, 2]

    @test ps.pb0.acols[1].nzval == T[4, 1, 1, -1, 1]
    @test ps.pb0.acols[2].nzval == T[4, 1, 1]
    @test ps.pb0.acols[3].nzval == T[4, 1, 3, -1]
    @test ps.pb0.acols[4].nzval == T[2, 1, 1]
    @test ps.pb0.acols[5].nzval == T[7]
    @test ps.pb0.acols[6].nzval == T[1]

    return nothing
end


function cg_strengthening_test2(T::Type)
    # We test CG_strengthening! function on the following MIP
    #=
    min     x + y + z
    s.t.    x + y         ⩾ 3.5
            x - y         ⩾ -2.9
                y + z     ⩽ 5.4
                y - z     ⩽ 7.25
            4 ⩾ x + y + z ⩾ 0

                x ⩾ 2
            5 ⩾ y ⩾ 1
            1 ⩾ z ⩾ -3

            x, y, z are integers
    =#

    C = T[1, 1, 1]
    lc = T[2, 1, -1]
    uc = T[Inf, 5, 1]
    lr = T[35//10, -29//10, -Inf, -Inf, 0]
    ur = T[Inf, Inf, 54//10, 725//100, 4]
    A = T[1 1 0
          1 -1 0
          0 1 1
          0 1 -1
          1 1 1]

    varTypes = [MOP.GENERAL_INTEGER, MOP.GENERAL_INTEGER,
                MOP.GENERAL_INTEGER]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(A), lr, ur, lc, uc,
        varTypes
    )
    ps = MOP.PresolveData(pb)

    MOP.cg_strengthening!(ps)

    @test ps.urow == T[Inf, Inf, 21, 14, 4]
    @test ps.lrow == T[4, -11, -Inf, -Inf, 0]

    @test ps.pb0.arows[1].nzval == T[1, 1]
    @test ps.pb0.arows[2].nzval == T[4, -4]
    @test ps.pb0.arows[3].nzval == T[4, 4]
    @test ps.pb0.arows[4].nzval == T[2, -2]

    @test ps.pb0.acols[1].nzval == T[1, 4, 1]
    @test ps.pb0.acols[2].nzval == T[1, -4, 4, 2, 1]
    @test ps.pb0.acols[3].nzval == T[4, -2, 1]
    return nothing
end


@testset "CG Strengthening Inequalities" begin
    for T in COEFF_TYPES
        @testset "$T" begin
            cg_strengthening_test1(T)
            cg_strengthening_test2(T)
        end
    end
end
