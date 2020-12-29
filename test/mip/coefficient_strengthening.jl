function positive_coef_strengthen_test(T::Type)
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
    lr = T[-1e30, -1e30]
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

    #MOP.coefficient_strenthening!(ps)

    #@test
    return ps
end

function negative_coef_strengthen_test(T::Type)
    Build the following model
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
   uc = T[2, 1, 2, 10^30]
   lr = T[-1e30, -1e30]
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

   #MOP.coefficient_strenthening!(ps)

   #@test
   return nothing
end
