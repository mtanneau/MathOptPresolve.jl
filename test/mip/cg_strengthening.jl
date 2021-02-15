function cg_strengthening_test1(T::Type)
    # We test CG_strengthening! function on the following MIP
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
    We will skip the the last 3 constraints because x₅ is continous
    coefficients of the forth constraint is not the same sign, and x₆
    is not bounded below
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
