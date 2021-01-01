function coef_strengthen_test1(T::Type)
    # Build the following model
    #=
        min     x + y + z + w
        s.t.     2 x - 5 y + 4 z + 3 w ⩽ 8
                         y +  z + 0.5w ⩽ 6
                 0 ⩽ x ⩽ 1
                 1 ⩽ y ⩽ 2.5
                 0 ⩽ z ⩽ 1
                 2 ⩽ w ⩽ 3
                 w is integer =#
    C = T[1, 1, 1, 1]
    lc = T[0, 1, 0, 2]
    uc = T[1, 5//2, 1, 3]
    lr = T[-Inf, -Inf]
    ur = T[8, 6]
    A = T[2 -5 4 3
          0 1 1 1]

    varTypes = [MOP.CONTINUOUS, MOP.CONTINUOUS,
                MOP.CONTINUOUS, MOP.GENERAL_INTEGER]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(A), lr, ur, lc, uc,
        varTypes
    )
    ps = MOP.PresolveData(pb)

    MOP.coefficient_strengthening!(ps)


    @test ps.urow == T[5, 9//2]
    @test ps.lrow == T[-Inf, -Inf]

    @test ps.pb0.arows[1].nzind == [1, 2, 3, 4]
    @test ps.pb0.arows[1].nzval == T[2, -5, 4, 2]
    @test ps.pb0.arows[2].nzind == [2, 3, 4]
    @test isapprox(ps.pb0.arows[2].nzval, T[1, 1, 1//2], atol = eps(T))

    @test ps.pb0.acols[1].nzind == [1]
    @test ps.pb0.acols[1].nzval == T[2]
    @test ps.pb0.acols[2].nzind == [1, 2]
    @test ps.pb0.acols[2].nzval == T[-5, 1]
    @test ps.pb0.acols[3].nzind == [1, 2]
    @test ps.pb0.acols[3].nzval == T[4, 1]
    @test ps.pb0.acols[4].nzind == [1, 2]
    @test ps.pb0.acols[4].nzval == T[2, 1//2]
    @test ps.status == MOP.NOT_INFERRED
    return nothing
end

function coef_strengthen_test2(T::Type)
   #=
       min     x + y + z + w
       s.t.     x + y + z - 2 w ⩽ 10
                2x + 3y + 4z + 5w ⩽ 20
                0 ⩽ x ⩽ 2
                0 ⩽ y ⩽ 1
                0 ⩽ z ⩽ 2
                -3 ⩽ w
                w is integer =#
   C = T[1, 1, 1, 1]
   lc = T[0, 0, 0, -3]
   uc = T[2, 1, 2, Inf]
   lr = T[-Inf, -Inf]
   ur = T[10, 20]
   A = T[1 1 1 -2
         2 3 4 5]

   varTypes = [MOP.CONTINUOUS, MOP.CONTINUOUS,
               MOP.CONTINUOUS, MOP.GENERAL_INTEGER]

   pb = MOP.ProblemData{T}()

   MOP.load_problem!(pb, "Test",
       true, C, zero(T),
       sparse(A), lr, ur, lc, uc,
       varTypes
   )
   ps = MOP.PresolveData(pb)

   MOP.coefficient_strengthening!(ps)

   @test ps.urow == T[7, 20]
   @test ps.lrow == T[-Inf, -Inf]

   @test ps.pb0.arows[1].nzind == [1, 2, 3, 4]
   @test ps.pb0.arows[1].nzval == T[1, 1, 1, -1]
   @test ps.pb0.arows[2].nzind == [1, 2, 3, 4]
   @test ps.pb0.arows[2].nzval == T[2, 3, 4, 5]

   @test ps.pb0.acols[1].nzind == [1, 2]
   @test ps.pb0.acols[1].nzval == T[1, 2]
   @test ps.pb0.acols[2].nzind == [1, 2]
   @test ps.pb0.acols[2].nzval == T[1, 3]
   @test ps.pb0.acols[3].nzind == [1, 2]
   @test ps.pb0.acols[3].nzval == T[1, 4]
   @test ps.pb0.acols[4].nzind == [1, 2]
   @test ps.pb0.acols[4].nzval == T[-1, 5]
   @test ps.status == MOP.NOT_INFERRED
   return nothing
end


function coef_strengthen_test3(T::Type)
    #=
        min     x₁ + x₂ + x₃ + x₄ + x₅
        s.t.     2 ⩽x₁ + -2x₂ + x₃ + x₄ + -x₅
                 x₁ + x₂ + x₄ + x₅ ⩽ 4.333
                 -4 ⩽ x₁ ⩽ 2
                 -1.5 ⩽ x₂ ⩽ 1
                 9 ⩽ x₃ ⩽ 15
                 0 ⩽ x₄ ⩽ 1
                 -1 ⩽ x₅ ⩽ 1
                 x₃ is integer, x₄ is binary  =#
    C = T[1, 1, 1, 1, 1]
    lc = T[-4, -3//2, 9, 0, -1]
    uc = T[2, 1, 15, 1, 1]
    lr = T[2, -Inf]
    ur = T[Inf, 13/3]
    A = T[1 -2 1 1 -1
          1 1 0 1 1]

    varTypes = [MOP.CONTINUOUS, MOP.CONTINUOUS,
                MOP.GENERAL_INTEGER, MOP.BINARY, MOP.CONTINUOUS]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(A), lr, ur, lc, uc,
        varTypes
    )
    ps = MOP.PresolveData(pb)

    MOP.coefficient_strengthening!(ps)


    @test ps.urow == T[Inf, 4]
    @test ps.lrow == T[-7, -Inf]

    @test ps.pb0.arows[1].nzind == [1, 2, 5]
    @test ps.pb0.arows[1].nzval == T[1, -2, -1]
    @test ps.pb0.arows[2].nzind == [1, 2, 4, 5]
    @test isapprox(ps.pb0.arows[2].nzval, T[1, 1, 2//3, 1], atol = 1e-12)


    @test ps.pb0.acols[1].nzind == [1, 2]
    @test ps.pb0.acols[1].nzval == T[1, 1]
    @test ps.pb0.acols[2].nzind == [1, 2]
    @test ps.pb0.acols[2].nzval == T[-2, 1]
    @test ps.pb0.acols[3].nzind == []
    @test ps.pb0.acols[3].nzval == T[]
    @test ps.pb0.acols[4].nzind == [2]
    @test isapprox(ps.pb0.acols[4].nzval, T[2//3], atol = 1e-12)

    @test ps.status == MOP.NOT_INFERRED
    return nothing
end

@testset "Coefficient Strengthening" begin
    for T in COEFF_TYPES
        @testset "$T" begin coef_strengthen_test1(T) end
        @testset "$T" begin coef_strengthen_test2(T) end
        @testset "$T" begin coef_strengthen_test3(T) end
    end
end
