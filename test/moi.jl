function _build_model(T::Type)
    # min  x[1]
    # s.t. x[1] >= -1
    #      -1 <= x[2] <= 3
    #      x[3] <= 1.2
    # .    x[3] in Z
    #      2.0 x[2] == 2.5
    model = MOIU.Model{T}()
    n = 3
    vis = MOI.add_variables(model, n)
    x = [MOI.SingleVariable(vis[i]) for i = 1:n]
    MOI.add_constraint(model, x[3], MOI.Integer())
    MOI.add_constraint(model, x[1], MOI.GreaterThan{T}(-1.0))
    MOI.add_constraint(model, x[2], MOI.Interval{T}(-1.0, 3.0))
    MOI.add_constraint(model, x[3], MOI.LessThan{T}(1.2))
    MOI.add_constraint(model, T(2.0) * x[2], MOI.EqualTo{T}(2.5))
    # TODO: Test when objective is not set
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), x[1])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return model
end

@testset "presolve!" begin
    for T in COEFF_TYPES
        src = _build_model(T)
        dest = MOIU.Model{T}()
        MOP.presolve!(dest, src, T)
        @test MOI.get(dest, MOI.NumberOfVariables()) == 0
        @test MOI.get(dest, MOI.ListOfConstraints()) == []
        @test MOI.get(dest, MOI.ObjectiveSense()) == MOI.MIN_SENSE
        @test MOI.get(dest, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}()) ==
              MOI.ScalarAffineFunction{T}([], -1.0)
    end
end
