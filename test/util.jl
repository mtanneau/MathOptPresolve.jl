function sync_columns_to_rows_test(T::Type)
    C = T[1, 0]
    lc = T[1, 0]
    uc = T[1, 0]
    lr = T[1]
    ur = T[0]
    A = T[1 2]

    varTypes = [MOP.GENERAL_INTEGER, MOP.CONTINUOUS]

    pb = MOP.ProblemData{T}()

    MOP.load_problem!(pb, "Test",
        true, C, zero(T),
        sparse(A), lr, ur, lc, uc,
        varTypes
    )
    ps = MOP.PresolveData(pb)

    # change ps.pb0.arows data
    ps.pb0.arows[1].nzval = [-1, 3]

    MOP.sync_columns_to_rows!(ps)

    @test ps.pb0.acols[1].nzval == T[-1]
    @test ps.pb0.acols[2].nzval == T[3]
    @test ps.pb0.acols[1].nzind == [1]
    @test ps.pb0.acols[2].nzind == [1]
    return nothing
end

@testset "Sync Columns to Rows" begin
    for T in COEFF_TYPES
        @testset "$T" begin
            sync_columns_to_rows_test(T)
        end
    end
end
