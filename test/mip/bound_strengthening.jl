function bound_strengthening_tests(T::Type)
end

@testset "Bound Strengthening" begin
    for T in COEFF_TYPES
        @testset "$T" begin bound_strengthening_tests(T) end
    end
end
