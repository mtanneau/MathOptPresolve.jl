var documenterSearchIndex = {"docs":
[{"location":"manual/reductions/","page":"Reductions","title":"Reductions","text":"CurrentModule = MathOptPresolve","category":"page"},{"location":"manual/reductions/#Primal-reductions","page":"Reductions","title":"Primal reductions","text":"","category":"section"},{"location":"manual/reductions/","page":"Reductions","title":"Reductions","text":"RemoveEmptyRow\nRemoveEmptyRows\nFixSingleVariable","category":"page"},{"location":"manual/reductions/#MathOptPresolve.RemoveEmptyRow","page":"Reductions","title":"MathOptPresolve.RemoveEmptyRow","text":"RemoveEmptyRow <: AbstractRule\n\nEliminate a single row if it is empty, i.e. all coefficients are zero.\n\nIn particular, consider row i given as\n\nl^r_i  sum_j a_i j x_j  u^r_i\n\nWe eliminate row i if a_ij = 0 for each variable index j.\n\nPresolve\n\nIf row i is empty, the corresponding constraint is completely removed from the problem.\n\nInfeasibility checks\n\nFor a prescribed tolerance ϵ  0, the problem is infeasible if either\n\nrow i is empty and its lower bound is positive, i.e., l^r_i  ϵ, or\nrow i is empty and its upper bound is negative, i.e., u^r_i  -ϵ.\n\nThe corresponding Farkas infeasibilities certificates are\n\ny^l = e_i, y^u = 0, z^l = z^u = 0,\ny^u = e_i, y^l = 0, z^l = z^u = 0,\n\nwhere e_i is the i-th vector of the canonical basis.\n\nPostsolve\n\nIf the problem is feasible, y^l_i = y^u_i = 0.\n\nIf the problem is infeasible, see the Farkas rays above.\n\nMisc\n\nThis is a primal reduction.\nThis reduction is non-destructive.\nThis reduction does not create any fill-in.\n\n\n\n\n\n","category":"type"},{"location":"manual/reductions/#MathOptPresolve.RemoveEmptyRows","page":"Reductions","title":"MathOptPresolve.RemoveEmptyRows","text":"RemoveEmptyRows <: AbstractRule\n\nRemove all empty rows in the problem.\n\nSee RemoveEmptyRow.\n\n\n\n\n\n","category":"type"},{"location":"manual/reductions/#MathOptPresolve.FixSingleVariable","page":"Reductions","title":"MathOptPresolve.FixSingleVariable","text":"FixSingleVariable <: AbstractRule\n\nEliminate a single variable x_j if its lower and upper bounds are equal.\n\nVariable x_j with lower and upper bounds l_j u_j is fixed to its lower bound if\n\n l_j - u_j  ϵ\n\nwhere ϵ is a prescribed tolerance.\n\nPresolve\n\nForall rows i,\n\nl^r_i  a_i jx_j + sum_k neq j a_i k x_k  u^r_i\n\nis transformed into\n\nl^r_i - a_i j l_j  sum_k neq j a_i k x_k  u^r_i - a_i j l_j\n\nThe objective constant c_0 is updated to c_0 + c_j l_j.\n\nInfeasibility checks\n\nThe problem is infeasible if a variable is fixed to a value outside its domain.\n\nGiven an integer feasibility tolerance ϵ_int  0, infeasibility is declared in the following cases:\n\nVariable x_j is integer (or implied integer) and l_j is fractional, i.e.,   min(l_j -  l_j   l_j  - l_j) geq ϵ_int.\nVariable x_j is binary, and   l_j notin -ϵ_int ϵ_int cup 1-ϵ_int 1+ϵ_int.\n\nPostsolve\n\nFor dual variables, first compute\n\nz_j = c_j - sum_i a_i j (y_i^l - y_i^u)\n\nthen recover z_j^l = z_j^+ and z_j^u = z_j^-.\n\nMisc\n\nThis is a primal reduction.\nThis reduction is non-destructive.\nThis reduction does not create any fill-in.\n\n\n\n\n\n","category":"type"},{"location":"reference/reference/","page":"Reference","title":"Reference","text":"Modules = [MathOptPresolve]\nOrder   = [:function]","category":"page"},{"location":"reference/reference/#MathOptPresolve.apply!","page":"Reference","title":"MathOptPresolve.apply!","text":"apply!(ps, rule, config) -> Nothing\n\nApply rule rule.\n\nArguments\n\nps::PresolveData{T}\nrule::AbstractRule\nconfig::PresolveOptions{T}\n\n\n\n\n\n","category":"function"},{"location":"reference/reference/#MathOptPresolve.bounds_consistency_checks!-Union{Tuple{MathOptPresolve.PresolveData{T}}, Tuple{T}} where T","page":"Reference","title":"MathOptPresolve.bounds_consistency_checks!","text":"bounds_consistency_checks(ps)\n\nCheck that all primal & dual bounds are consistent.\n\nTODO: If not, declare primal/dual infeasibility and extract ray.\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.load_problem!-Union{Tuple{T}, Tuple{MathOptPresolve.ProblemData{T},String,Bool,Array{T,1},T,SparseArrays.SparseMatrixCSC,Array{T,1},Array{T,1},Array{T,1},Array{T,1}}, Tuple{MathOptPresolve.ProblemData{T},String,Bool,Array{T,1},T,SparseArrays.SparseMatrixCSC,Array{T,1},Array{T,1},Array{T,1},Array{T,1},Union{Nothing, Array{MathOptPresolve.VariableType,1}}}} where T","page":"Reference","title":"MathOptPresolve.load_problem!","text":"load_problem!(pb, ...)\n\nLoad entire problem.\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.postsolve!-Union{Tuple{T}, Tuple{MathOptPresolve.Solution{T},MathOptPresolve.Solution{T},MathOptPresolve.PresolveData{T}}} where T","page":"Reference","title":"MathOptPresolve.postsolve!","text":"postsolve!\n\nPerform post-solve.\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.presolve!-Union{Tuple{MathOptPresolve.PresolveData{T}}, Tuple{T}} where T","page":"Reference","title":"MathOptPresolve.presolve!","text":"presolve(pb::ProblemData)\n\nPerform pre-solve.\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.remove_empty_columns!-Union{Tuple{MathOptPresolve.PresolveData{T}}, Tuple{T}} where T","page":"Reference","title":"MathOptPresolve.remove_empty_columns!","text":"remove_empty_columns!(ps::PresolveData)\n\nRemove all empty columns.\n\nCalled once at the beginning of the presolve procedure. If an empty column is created later, it is removed on the spot.\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.remove_forcing_rows!-Tuple{MathOptPresolve.PresolveData}","page":"Reference","title":"MathOptPresolve.remove_forcing_rows!","text":"remove_forcing_rows!\n\nRemove forcing and dominated row\n\n\n\n\n\n","category":"method"},{"location":"reference/reference/#MathOptPresolve.remove_free_column_singletons!-Tuple{MathOptPresolve.PresolveData}","page":"Reference","title":"MathOptPresolve.remove_free_column_singletons!","text":"remove_free_column_singletons!(ps)\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MathOptPresolve","category":"page"},{"location":"#MathOptPresolve","page":"Home","title":"MathOptPresolve","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"warning: Warning\nThe documention is under construction, and should not be considered complete at this time.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
