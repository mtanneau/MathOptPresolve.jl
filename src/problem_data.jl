using SparseArrays

mutable struct RowOrCol{T}
    nzind::Vector{Int}
    nzval::Vector{T}
end

const Row = RowOrCol
const Col = RowOrCol

"""
    ProblemData{T}

Data structure for storing problem data in precision `T`.

The LP is represented in canonical form

```math
\\begin{array}{rl}
    \\displaystyle \\min_{x} \\ \\ \\ & c^{T} x + c_{0} \\\\
    s.t. \\ \\ \\ & l_{r} \\leq A x \\leq u_{r} \\\\
    & l_{c} \\leq x \\leq u_{c}
\\end{array}
```
"""
mutable struct ProblemData{T}

    name::String

    # Dimensions
    ncon::Int  # Number of rows
    nvar::Int  # Number of columns (i.e. variables)

    # Objective
    # TODO: objective sense
    objsense::Bool  # true is min, false is max
    obj::Vector{T}
    obj0::T  # Constant objective offset

    # Constraint matrix
    # We store both rows and columns. It is redundant but simplifies access.
    # TODO: put this in its own data structure? (would allow more flexibility in modelling)
    arows::Vector{Row{T}}
    acols::Vector{Col{T}}

    # TODO: Data structures for QP
    # qrows
    # qcols
    
    # Bounds
    lcon::Vector{T}
    ucon::Vector{T}
    lvar::Vector{T}
    uvar::Vector{T}

    # Names
    con_names::Vector{String}
    var_names::Vector{String}

    # Only allow empty problems to be instantiated for now
    ProblemData{T}(pbname::String="") where {T} = new{T}(
        pbname, 0, 0,
        true, T[], zero(T),
        Row{T}[], Col{T}[],
        T[], T[], T[], T[],
        String[], String[]
    )
end

import Base.empty!

function Base.empty!(pb::ProblemData{T}) where {T}

    pb.name = ""

    pb.ncon = 0
    pb.nvar = 0

    pb.objsense = true
    pb.obj = T[]
    pb.obj0 = zero(T)

    pb.arows = Row{T}[]
    pb.acols = Col{T}[]

    pb.lcon = T[]
    pb.ucon = T[]
    pb.lvar = T[]
    pb.uvar = T[]

    pb.con_names = String[]
    pb.var_names = String[]

    return pb
end

"""
    load_problem!(pb, )

Load entire problem.
"""
function load_problem!(pb::ProblemData{T},
    name::String,
    objsense::Bool, obj::Vector{T}, obj0::T,
    A::SparseMatrixCSC,
    lcon::Vector{T}, ucon::Vector{T},
    lvar::Vector{T}, uvar::Vector{T},
    con_names::Vector{String}, var_names::Vector{String}
) where {T}
    empty!(pb)

    # Sanity checks
    ncon, nvar = size(A)
    ncon == length(lcon) || error("")
    ncon == length(ucon) || error("")
    ncon == length(con_names) || error("")
    nvar == length(obj)
    isfinite(obj0) || error("Objective offset $obj0 is not finite")
    nvar == length(lvar) || error("")
    nvar == length(uvar) || error("")

    # Copy data
    pb.name = name
    pb.ncon = ncon
    pb.nvar = nvar
    pb.objsense = objsense
    pb.obj = copy(obj)
    pb.obj0 = obj0
    pb.lcon = copy(lcon)
    pb.ucon = copy(ucon)
    pb.lvar = copy(lvar)
    pb.uvar = copy(uvar)
    pb.con_names = copy(con_names)
    pb.var_names = copy(var_names)

    # Load coefficients
    pb.acols = Vector{Col{T}}(undef, nvar)
    pb.arows = Vector{Row{T}}(undef, ncon)
    for j in 1:nvar
        col = A[:, j]
        pb.acols[j] = Col{T}(col.nzind, col.nzval)
    end

    At = sparse(A')
    for i in 1:ncon
        row = At[:, i]
        pb.arows[i] = Row{T}(row.nzind, row.nzval)
    end

    return pb
end
