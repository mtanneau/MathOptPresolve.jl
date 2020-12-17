function ensure_integer_bounds!(ps::PresolveData, j::Int)
    # Column was already removed.
    ps.colflag[j] || return nothing

    if (ps.var_types[j] == BINARY)
        ps.lcol[j] = T(max(0, ceil(ps.lcol[j])))
         ps.ucol[j] = T(min(1, floor(ps.ucol[j])))
    elseif (ps.var_types[j] == GENERAL_INTEGER)
        ps.lcol[j], ps.ucol[j] = T(ceil(ps.lcol[j])), T(floor(ps.ucol[j]))
    end
    return nothing
end
