# Positive and negative part of a number
pos_part(x::T) where {T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where {T} = x >= zero(T) ? zero(T) : -x

# Allow some numerical errors to avoid the situation
# like ceiling 1e-10 to 1.
function approx_ceil(val::T, ϵ_int::T)::T where{T}
    if (val - floor(val)) < ϵ_int
        return T(floor(val))
    end
    return T(ceil(val))
end

# Allow some numerical errors to avoid the situation
# like flooring 1 - 1e-10 to 0.
function approx_floor(val::T, ϵ_int::T)::T where{T}
    if (ceil(val) - val) < ϵ_int
        return T(ceil(val))
    end
    return T(floor(val))
end
