# Bipartite matching, using the Ford–Fulkerson algorithm (with DFS). (Algorithms
# with better asymptotic running time exist, but they may not be more efficient
# in practice, and they tend to be quite a bit more complex.) `X` should be a
# matrix where `X[i, g] ≠ 0` indicates an edge between `i` and `j`. I.e., this
# may be a boolean matrix, or any other numeric matrix where a matching along
# nonzero values (e.g., using the valuation matrix directly in the case of MNW)
# is desired. The function returns an iterator over pairs `(i => g)`,
# indicating that `i` is matched with `g`.
function bipartite_matching(X)

    n, m = size(X)

    mate, seen = zeros(Int, m), falses(m)

    function match!(i)
        for g = 1:m
            (iszero(X[i, g]) || seen[g]) && continue
            seen[g] = true
            if mate[g] == 0 || match!(mate[g])
                mate[g] = i
                return true
            end
        end
        return false
    end

    for i = 1:n
        seen .= false
        match!(i)
    end

    return ((i => g) for (g, i) in pairs(mate) if i ≠ 0)

end


bipartite_matching(X::Profile) = bipartite_matching(matrix(X))


# Bare-bones lexicographic optimization. The model is a JuMP model, and the
# objectives are given as a sequence of `(sense, func)` pairs, acting like the
# parameters `sense` and `func` in `JuMP.@objective`. The (absolute) tolerance ϵ
# is the same as in `solve_mip`, and is used when adding constraints at each
# step.
#
# Note: No cleanup is performed, in order for `JuMP.value` to function properly
# on the model afterward (i.e., to avoid a `JuMP.OptimizeNotCalled` exception).
# This means that the n–1 constraints added, to lock in all but the last
# objective, will all still be there after this function has been run.
# In order to permit "manual" cleanup, a vector of the added constraints is
# returned.
function lex_optimize!(model, objectives; ϵ=1e-5)

    constraints = Any[]

    nobj = length(objectives)

    for (i, (sense, func)) in enumerate(objectives)

        @objective(model, sense, func)
        optimize!(model)
        status = termination_status(model)

        i == nobj && break

        status in conf.MIP_SUCCESS || break

        # Lock in this objective, before going to the next one:
        z = objective_value(model)
        if sense == MOI.MIN_SENSE
            con = @constraint(model, func <= z + ϵ)
        else
            con = @constraint(model, func >= z - ϵ)
        end

        push!(constraints, con)

    end

    return constraints

end


# Takes in a model and a set of functions (JuMP affine expressions) that are to
# be turned into a leximax objective (i.e., first minimize the maximum, then the
# second largest value, etc.). Then constructs a set of objective functions
# to use instead, which should be minimized lexicographically (first minimize
# the first, then the second, etc.). This is a direct implementation of the
# scheme described by Theorem 1 of Ogryczak and Śliwiński in their paper ["On
# Direct Methods for Lexicographic Min-Max
# Optimization"](https://doi.org/10.1007/11751595_85) (2006). The main point of
# this construction is that it works on nonconvex (e.g., discrete) problems.
#
# Note: They conclude that their second approach is superior, but that is based
# on a limited set of objective values, which may be too strong a requirement in
# our case.
#
# The function also modifies the model in-place, adding some variables and
# constraints. For readability, named variables `t` and `d` are used, as in the
# paper. If these collide with variables used in the main program, simply
# specify other names using the `tname` and `dname` keyword arguments.
#
# In order to get leximin instead, simply negate the functions. (The resulting
# objectives should then still be minimized.)
function leximax_os06_1!(model, funcs; tname="t", dname="d")

    m = length(funcs)

    @variable(model, t[1:m],           base_name=tname)
    @variable(model, d[1:m, 1:m] >= 0, base_name=dname)

    for j = 1:m, k = 1:m
        @constraint(model, t[k] + d[k, j] >= funcs[j])
    end

    return [i * t[i] + sum(d[i,j] for j=1:m) for i = 1:m]

end
