# Positive and negative part of a number
pos_part(x::T) where {T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where {T} = x >= zero(T) ? zero(T) : -x


function calc_upper_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        if (val > 0)
            isfinite(ucol[ind]) || return Inf
            bound += val * ucol[ind]
        else
            isfinite(lcol[ind]) || return Inf
            bound += val * lcol[ind]
        end
    end
    return bound
end

function calc_lower_bound_except_one(row::Row{T}, lcol::Vector{T},
                            ucol::Vector{T}, j::Int) where {T}
    bound = T(0)
    for i in 1:length(row.nzind)
        ind, val = row.nzind[i], row.nzval[i]
        (ind == j) && continue
        if (val > 0)
            isfinite(lcol[ind]) || return -Inf
            bound += val * lcol[ind]
        else
            isfinite(ucol[ind]) || return -Inf
            bound += val * ucol[ind]
        end
    end
    return bound
end

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
