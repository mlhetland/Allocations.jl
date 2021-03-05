# Algorithms that share some allocation-specific infrastructure based on
# mixed-integer programming using JuMP.


# A shared context for the MIP pipeline. The JuMP variable representing the
# allocation (as an n-by-m matrix) is kept in alloc_var. The res field should be
# a named tuple of any extra data to be splatted in at the end of the result.
mutable struct MIPContext{V <: Additive, M, A, S}
    valuation::V
    alloc_var::A
    model::M
    solver::S
    alloc::Union{Allocation, Nothing}
    res::NamedTuple
end
MIPContext(v, a, m, s) = MIPContext(v, a, m, s, nothing, (;))


# Internal MIP-building functions

# The init_mip and solve_mip functions are used to set up a MIP based on a
# valuation, and then to solve the MIP and produce an Allocation. After this,
# a final stage (such as mnw_result) should be used, to wrap the allocation
# and any other result variables (such as the objective value) in a named
# tuple.
#
# The "enforce" functions add constraints to the MIP, while the "achieve"
# functions also modify the objective function. At most one "achieve" function
# should be used in a single pipeline; otherwise, the second one would
# overwrite the objective of the first one.


# Set up a MIPContext with based on a valuation. Initializes a JuMP model with
# constraints ensuring a valid allocation.
function init_mip(V::Additive, solver)

    model = Model(solver)

    N = agents(V)
    M = items(V)

    @variable(model, A[N, M], binary=true)

    for g in M
        # Each item allocated to at exactly one agent
        @constraint(model, sum(A[i, g] for i in N) == 1)
    end

    return MIPContext(V, A, model, solver)

end


# Shortcut, so you can use a matrix instead of an Additive valuation,
# including in the alloc_... functions that rely on init_mip.
init_mip(V::Matrix, solver) = init_mip(Additive(V), solver)


# Solves the MIP model and constructs the actual Allocation object in the
# MIPContext.
function solve_mip(ctx)

    optimize!(ctx.model)

    @assert termination_status(ctx.model) in conf.MIP_SUCCESS

    V = ctx.valuation
    ctx.alloc = Allocation(na(V), ni(V))

    for i in agents(V), g in items(V)
        JuMP.value(ctx.alloc_var[i, g]) ≈ 1.0 && give!(ctx.alloc, i, g)
    end

    return ctx

end


## Objective-modifying pipeline steps (achieve_...)


# Set up objective and constraints to make sure the JuMP model produces an MNW
# allocation, using a slightly adapted version of the approach of Caragiannis
# et al. (https://doi.org/10.1145/3355902).
achieve_mnw(mnw_warn) = function(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    @assert isintegral(V)
    @assert isnonnegative(V)

    positive = [value(V, i, j) > 0 for i in agents(V), j in items(V)]

    matching = bipartite_matching(positive, ctx.solver)

    N = findall(vec(any(matching, dims=2)))

    M = items(V)

    v_max = Float64(maximum(value(V, i, M) for i in N))

    mnw_prec = true

    for (k, name) in [1 => "PO", 2 => "EF1", length(N) => "MNW"]
        if log(v_max^k) - log(v_max^k - 1) == 0.0
            if name == "MNW"
                mnw_prec = false
                !mnw_warn && continue
            end
            @warn("Precision insufficient to guarantee $name")
        end
    end

    ctx.res = (ctx.res..., mnw_prec = mnw_prec)

    @variable(model, W[N])

    @objective(model, Max, sum(W))

    A = ctx.alloc_var

    for i in N, k = 1:2:v_max
        @constraint(model, W[i] <=
                log(k) + (log(k + 1) - log(k)) *
                (sum(A[i, g] * value(V, i, g) for g = M) - k))
    end

    for i in N
        @constraint(model, sum(A[i, g] * value(V, i, g) for g in M) >= 1)
    end

    return ctx

end


# Set up objective and constraints to make sure the JuMP model produces an
# egalitarian/maximin allocation.
achieve_mm(cutoff=nothing) = function(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    @variable(model, v_min)

    for i in N
        @constraint(model, v_min <= sum(A[i, g] * value(V, i, g) for g in M))
    end

    if cutoff ≢ nothing
        @constraint(model, v_min <= cutoff)
    end

    @objective(model, Max, v_min)

    return ctx

end


# Set up objective and constraints to make sure the JuMP model produces a
# minimum ordered weighted average, by utility rank, using wt as the weight
# function (cf., Lesca & Perny, 2010)
achieve_mgg(wt) = function(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    N, n, M = agents(V), na(V), items(V)

    @variable(model, r[N])
    @variable(model, b[N, N] >= 0)

    omega = [i == n ? wt(i, n) : wt(i, n) - wt(i + 1, n) for i in N]

    @objective(model, Max,
        sum(omega[j] * (j * r[j] - sum(b[i, j] for i in N)) for j in N))

    for i in N, j in N
        @constraint(model,
            r[j] - b[i, j] <= sum(A[i, g] * value(V, i, g) for g in M))
    end

    return ctx

end


## Objective-preserving pipeline steps (enforce_...)


# Enforce no constraints on the JuMP model.
enforce(C::Nothing) = identity


# Enforce cardinality constraints on the JuMP model.
enforce(C::Counts) = function(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    for S in C, i in agents(V)
        @constraint(model, sum(A[i, g] for g in S) <= S.threshold)
    end

    return ctx

end


# Enforce item conflict constraints on the JuMP model.
enforce(C::Conflicts) = function(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    G = graph(C)

    @assert nv(G) == ni(V)

    for e in edges(G)
        g, h = src(e), dst(e)
        for i in agents(V)
            @constraint(model, A[i, g] + A[i, h] <= 1)
        end
    end

    return ctx

end


# Enforce envy-freeness up to a single object (EF1) on the JuMP model.
function enforce_ef1(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    # D[i, j, g]: Item g dropped for i to not envy j
    @variable(model, D[N, N, M], binary=true)

    for i in N, j in N

        # Drop at most one item for envy-freeness
        @constraint(model, sum(D[i, j, g] for g in M) <= 1)

        for g in M
            # Drop only items allocated to j to make i not envy j
            @constraint(model, D[i, j, g] <= A[j, g])
        end

    end

    for i in N, j in N

        i == j && continue

        # Agent i's value for her own bundle
        vii = sum(value(V, i, g) * A[i, g] for g in M)

        # Agent i's value for j's bundle, without dropped item
        vij1 = sum(value(V, i, g) * (A[j, g] - D[i, j, g]) for g in M)

        # No envy, once an item is (possibly) dropped:
        @constraint(model, vii >= vij1)

    end

    return ctx

end


# Enforce envy-freeness up to any object (EFX) on the JuMP model.
function enforce_efx(ctx)

    @warn "Buggy -- currently equivalent to EF"

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    for i in N, j in N

        i == j && continue

        # Agent i's value for her own bundle
        vii = sum(value(V, i, g) * A[i, g] for g in M)

        # XXX This doesn't work, as it'll also include constraints where we skip
        # objects outside j's bundle:

        # For each possibly dropped item d ...
        for d in M

            # Agent i's value for j's bundle, without d
            vijx = sum(value(V, i, g) * A[j, g] for g in M if g ≠ d)

            # No envy, once d is dropped
            @constraint(model, vii >= vijx)

        end

    end

    return ctx

end


# Actual allocation methods


# Generic function to extract the allocation at the end of the pipeline. Many
# allocation methods may want their own `..._result` functions, adding other
# fields to the named tuple being returned.
alloc_result(ctx) = (alloc=ctx.alloc, model=ctx.model, ctx.res...)


"""
    alloc_ef1(V, C; solver=conf.MIP_SOLVER)

Create an `Allocation` that is envy-free up to one item (EF1), based on the
valuation `V`, possibly subject to the constraints given by the `Constraint`
object `C`. The solution is found using a straightforward mixed-integer program,
and is most suitable for constraints where no specialized algorithm exists. For
example, without constraints, a straightforward round robin picking sequence
yields EF1, and a similar strategy works for cardinality constraints. (It is
still possible to use this function without constraints, by explicitly supplying
`nothing` for the constraint argument `C`.) The return value is a named tuple
with the fields `alloc` (the `Allocation`) and `model` (the JuMP model used in
the computation).

Note that for some constraints, there may not *be* an EF1 allocation, in which
case the function will fail with an exception.
"""
function alloc_ef1(V, C; solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    enforce_ef1 |>
    enforce(C) |>
    solve_mip |>
    alloc_result

end


"""
    alloc_efx(V[, C]; solver=conf.MIP_SOLVER)

Create an `Allocation` that is envy-free up to any item (EFX), based on the
valuation `V`, possibly subject to the constraints given by the `Constraint`
object `C`. The solution is found using a straightforward mixed-integer program.
The return value is a named tuple with the fields `alloc` (the `Allocation`) and
`model` (the JuMP model used in the computation).

Note that while some constraints may prevent an exact EFX allocation, it is
currently (Mar 2021) an open question whether EFX always exists in the
unconstrained case (see, e.g., [*Improving EFX Guarantees through Rainbow Cycle
Number*](https://arxiv.org/abs/2103.01628) by Chaudhury et al.).
"""
function alloc_efx(V, C=nothing; solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    enforce_efx |>
    enforce(C) |>
    solve_mip |>
    alloc_result

end


# Extract the allocation and the MNW value (excluding agents with a utility of
# zero) at the end of the pipeline. Strictly speaking, we needn't include the
# mnw field here, as it's a separate calculations; however, it's convenient if
# you supply a matrix as the argument for alloc_mnw function, since we have
# access to the valuation here. (Also, the cost of nash_welfare will generally
# be very low, compared to the actual MIP solving.)
function mnw_result(ctx)
    V, A = ctx.valuation, ctx.alloc
    return (alloc=A, model=ctx.model, mnw=nash_welfare(V, A), ctx.res...)
end


"""
    alloc_mnw(V[, C]; mnw_warn=true, solver=conf.MIP_SOLVER)

Create an `Allocation` attaining maximum Nash welfare (MNW), based on the
valuation `V`, possibly subject to the constraints given by the `Constraint`
object `C`. The solution is found using the approach of Caragiannis et al. in
their 2019 paper [The Unreasonable Fairness of Maximum Nash
Welfare](https://doi.org/10.1145/3355902), with two minor modifications:

1. Rather than hard-coding a maximum valuation (arising from the assumption
   that the values of each agent sum to 1000), this maximum is extracted from
   `V`; and

2. Extra constraints are permitted (through the object `C`), possibly lowering
   the attainable MNW.

Because of how the integer program is constructed, it is sensitive to
precision effects, where a high number of agents, can make it impossible to
guarantee Pareto optimalty (PO), EF1 or MNW respectively. If the precision is
too low, the appropriate warning will be issued, but the computation is not
halted. It may be useful to find a solution that is guaranteed to satisfy PO and
EF1, even if it may not be exactly MNW. For such cases, the `mnw_warn` keyword
argument may be set to `false`, to suppress the MNW warning.

The return value is a named tuple with fields `alloc` (the `Allocation`),
`mnw` (the achieved Nash welfare for the agents with nonzero utility),
`mnw_prec` (whether or not there was enough precision to find MNW) and `model`
(the JuMP model used in the computation).
"""
function alloc_mnw(V, C=nothing; mnw_warn=true, solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mnw(mnw_warn) |>
    enforce(C) |>
    solve_mip |>
    mnw_result

end


"""
    alloc_mnw_ef1(V, C; mnw_warn=true, solver=conf.MIP_SOLVER)

Equivalent to `alloc_mnw`, except that EF1 is enforced. Without any added
constraints, MNW implies EF1, so this function is not needed in that case.
Therefore the argument `C` is not optional.
"""
function alloc_mnw_ef1(V, C; mnw_warn=true, solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mnw(mnw_warn) |>
    enforce_ef1 |>
    enforce(C) |>
    solve_mip |>
    mnw_result

end


# Extract the allocation and the maximin value at the end of the pipeline.
function mm_result(ctx)
    return (alloc=ctx.alloc,
            model=ctx.model,
            mm=objective_value(ctx.model),
            ctx.res...)
end


"""
    alloc_mm(V[, C]; cutoff=nothing, solver=conf.MIP_SOLVER)

Create an egalitarion or maximin `Allocation`, i.e., one where the minimum
bundle value is maximized. The `cutoff`, if any, is a level at which we are
satisfied, i.e., any allocation where all agents attain this value is
acceptable. The return value is a named tuple with fields `alloc` (the
`Allocation`), `mm` (the lowest agent utility) and `model` (the JuMP model
used in the computation).
"""
function alloc_mm(V, C=nothing; cutoff=nothing, solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mm(cutoff) |>
    enforce(C) |>
    solve_mip |>
    mm_result

end


"""
    mms(V, i[, C]; solver=conf.MIP_SOLVER)

Determine the maximin share of agent `i`, i.e., the bundle value she is
guaranteed to attain if she partitions the items and the other agents choose
their bundles. Useful, e.g., as a point of reference when determining the
empirical approximation ratios of approximate MMS allocation algorithms. Also
used as a subroutine in `alloc_mms`. The return value is a named tuple with the
fields `mms` (the maximin share of agent `i`) and `model` (the JuMP model used
in the computation).
"""
function mms(V::Additive, i, C=nothing; solver=conf.MIP_SOLVER)

    # Let all agents be clones of agent i
    Vi = Additive([value(V, i, g) for _ in agents(V), g in items(V)])

    # maximin in this scenario is MMS for agent i
    res = alloc_mm(Vi, C, solver=solver)

    return (mms=res.mm, model=res.model)

end


"""
    mms_alpha(V, A, mmss)

Utility function to find the fraction of the maximin share guarantee attained by
the allocation `A`, under the valuation `V`, where `mmss[i]` is the MMS of agent
`i`. This makes it possible, for example, to use the `mmss` field from the
result of `alloc_mms` to find the MMS approximation provided by an allocation
constructed by other means. For example:

    mmss = alloc_mms(V).mmss
    A = alloc_rand(V).alloc
    alpha = mms_alpha(V, A, mmss)

"""
function mms_alpha(V, A, mmss)
    vals = [value(V, i, A) for i in agents(A)]
    return minimum(vals ./ mmss)
end


"""
    alloc_mms(V[, C]; cutoff=false, solver=conf.MIP_SOLVER)

Find an MMS allocation, i.e., one that satisfies the *maximin share
guarantee*, where each agent gets a bundle it weakly prefers to its maximin
share (introduced by Budish, in his 2011 paper [The Combinatorial Assignment
Problem: Approximate Competitive Equilibrium from Equal
Incomes](https://doi.org/10.1086/664613)). The return value is a named tuple
with fields `alloc` (the `Allocation`), `mmss`, the individual MMS values for
the instance, `alpha`, the lowest fraction of MMS that any agent achieves
(is at least 1 exactly when the allocation is MMS), `model` (the JuMP model used
in computing `alpha`) and `mms_models` (the JuMP models used to compute the
individual maximin shares). If `cutoff` is set to `true`, this fraction is
capped at 1.
"""
function alloc_mms(V::Additive, C=nothing; cutoff=false, solver=conf.MIP_SOLVER)

    N, M = agents(V), items(V)

    X = zeros(na(V), ni(V))

    ress = [mms(V, i, C, solver=solver) for i in N]

    # individual MMS values -- also included in the result
    mmss = [res.mms for res in ress]

    mms_models = [res.model for res in ress]

    for i in N

        for g in M
            X[i, g] = value(V, i, g) / mmss[i]
        end

    end

    max_alpha = cutoff ? 1.0 : nothing

    # maximin with scaled values is as close to the MMS guarantee as possible
    res = alloc_mm(Additive(X), C, cutoff=max_alpha, solver=solver)

    return (alloc=res.alloc, model=res.model, mms_models=mms_models,
            alpha=res.mm, mmss=mmss)

end


alloc_mms(V::Matrix, C=nothing; cutoff=false, solver=conf.MIP_SOLVER) =
    alloc_mms(Additive(V), C, cutoff=cutoff, solver=solver)


# Extract the allocation at the end of the pipeline.
function mgg_result(ctx)
    return (alloc=ctx.alloc, model=ctx.model, ctx.res...)
end


"""
    wt_gini(n, i)

The (unnormalized) weights used in the ordered weighted average in the Gini
social-evaluation function, where the utility of the `i`th agent, ordered by
increasing utility, is given weight ``2(n - i) + 1``. (The normalized weights
yielding the original Gini social-evaluation function are divided by ``n^2``,
but this makes no difference to the optimization problem.)
"""
wt_gini(i, n) = 2(n - i) + 1


"""
    alloc_mgg(V[, C]; wt=wt_gini, solver=conf.MIP_SOLVER)

Minimizes the generalized Gini social-evaluation function, or other ordered
weighted averages of agent utilities, where the weight is based on utility
rank `i`, from the least happy (`1`) to the most happy (`n`), parameterized by
the function `wt(i, n)`. The method used is based on that of Lesca and Perny
(linear formulation ``\\Pi'_W``, without ``\\alpha``, ``\\alpha'``, ``\\beta``
and ``\\beta'``) in their paper 2010 paper [LP Solvable Models for Multiagent
Fair Allocation
Problems](https://pdfs.semanticscholar.org/c6ae/41213e5744ec0a1cca632155a42a46f8b2ad.pdf).
The return value is a named tuple with the fields `alloc` (the `Allocation`
that has been produced) and `model` (the JuMP model used in the computation).
"""
function alloc_mgg(V, C=nothing; wt=wt_gini, solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mgg(wt) |>
    enforce(C) |>
    solve_mip |>
    mgg_result

end
