function ensure_integer_bounds(ps::PresolveData, j::Int)
    # Column was already removed or
    ps.colflag[j] || (ps.var_types[j] != CONTINUOUS) || return nothing

    ps.lcol[j], ps.ucol[j] = ceil(ps.lcol[j]), floor(ps.ucol[j])
    return nothing
end
