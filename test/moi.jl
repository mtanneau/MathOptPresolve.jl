function _is_approx_integral(x::Real)
    x_c = ceil(Int, x)
    x_f = floor(Int, x)
    return abs(x - x_c) ≈ 0 || abs(x - x_f) ≈ 0
end

@testset "presolve!" begin
    for T in COEFF_TYPES
        @testset "solved" begin
            src = MOIU.Model{T}()
            n = 3
            x = MOI.add_variables(src, n)
            MOI.add_constraint(src, x[3], MOI.Integer())
            MOI.add_constraint(src, x[1], MOI.GreaterThan{T}(-1.0))
            MOI.add_constraint(src, x[2], MOI.Interval{T}(-1.0, 3.0))
            MOI.add_constraint(src, x[3], MOI.LessThan{T}(1.2))
            MOI.add_constraint(src, T(2.0) * x[2], MOI.EqualTo{T}(2.5))
            MOI.set(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), one(T) * x[1])
            MOI.set(src, MOI.ObjectiveSense(), MOI.MIN_SENSE)

            dest = MOIU.Model{T}()
            pr = @inferred MOP.presolve!(dest, src, T)

            @test get_status(pr) == MOI.OPTIMAL
            for x in (get_optimal_solution(pr), post_crush(pr, T[]))
                @test length(x) == 3
                @test x[1] ≈ T(-1.0)
                @test x[2] ≈ T(2.5 / 2.0)
                @test x[3] <= 1.2
                @test _is_approx_integral(x[3])
            end

            @test_throws ArgumentError post_crush(pr, T[1.0])
            @test_throws ErrorException get_unbounded_ray(pr)
            @test_throws ErrorException get_infeasibility_certificate(pr)
        end
        @testset "infeasible" begin
            src = MOIU.Model{T}()
            n = 5
            x = MOI.add_variables(src, n)
            MOI.add_constraint(src, x[1], MOI.Interval{T}(0.0, 1.0))
            MOI.add_constraint(src, T(3.0) * x[1] + T(3.0), MOI.LessThan{T}(0.0))

            dest = MOIU.Model{T}()
            pr = @inferred MOP.presolve!(dest, src, T)

            @test get_status(pr) == MOI.INFEASIBLE

            @test_throws ArgumentError post_crush(pr, T[])
            @test_throws ErrorException get_optimal_solution(pr)
            @test_throws ErrorException get_unbounded_ray(pr)
        end
        @testset "unbounded" begin
            src = MOIU.Model{T}()
            x = MOI.add_variable(src)
            MOI.add_constraint(src, x, MOI.GreaterThan{T}(2.0))
            MOI.set(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), one(T) * x)
            MOI.set(src, MOI.ObjectiveSense(), MOI.MAX_SENSE)

            dest = MOIU.Model{T}()
            pr = @inferred MOP.presolve!(dest, src, T)

            @test get_status(pr) == MOI.DUAL_INFEASIBLE
            @test_throws ErrorException get_optimal_solution(pr)
            @test get_unbounded_ray(pr) ≈ T[1.0]

            @test post_crush(pr, T[3.4], is_ray = true) ≈ T[3.4]
        end
        @testset "not inferred" begin
            src = MOIU.Model{T}()
            n = 3
            x = MOI.add_variables(src, n)
            MOI.add_constraint(src, T(1.0) * x[1] + T(1.0) * x[2], MOI.GreaterThan{T}(0.0))
            MOI.add_constraint(src, T(1.0) * x[1] - T(1.0) * x[2], MOI.GreaterThan{T}(0.0))
            MOI.add_constraint(src, T(2.3) * x[3], MOI.EqualTo{T}(2.3))
            MOI.set(
                src,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
                T(1.0) * x[2] + T(1.0) * x[3],
            )
            MOI.set(src, MOI.ObjectiveSense(), MOI.MIN_SENSE)

            dest = MOIU.Model{T}()
            pr = @inferred MOP.presolve!(dest, src, T)

            @test get_status(pr) == MOI.OPTIMIZE_NOT_CALLED
            @test_throws ErrorException get_optimal_solution(pr)
            @test_throws ErrorException get_unbounded_ray(pr)
        end
    end
end
