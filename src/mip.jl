# Algorithms that share some allocation-specific infrastructure based on
# mixed-integer programming using JuMP.


# A shared context for the MIP pipeline. The JuMP variable representing the
# allocation (as an n-by-m matrix) is kept in alloc_var. The res field should be
# a named tuple of any extra data to be splatted in at the end of the result.
# The objectives field is only used with lexicographic optimization, where
# solve_mip calls lex_optimize!.
mutable struct MIPContext{V <: Additive, M, A, S}
    profile::V
    alloc_var::A
    model::M
    solver::S
    objectives::Vector{Any}
    alloc::Union{Allocation, Nothing}
    callbacks::Vector{Function}
    res::NamedTuple
end
MIPContext(v, a, m, s) = MIPContext(v, a, m, s, [], nothing, Function[], (;))


na(ctx::MIPContext) = na(ctx.profile)
ni(ctx::MIPContext) = ni(ctx.profile)


# Internal MIP-building functions

# The init_mip and solve_mip functions are used to set up a MIP based on a
# valuation profile, and then to solve the MIP and produce an Allocation. After
# this, a final stage (such as mnw_result) should be used, to wrap the
# allocation and any other result variables (such as the objective value) in a
# named tuple.
#
# The "enforce" functions add constraints to the MIP, while the "achieve"
# functions also modify the objective function. At most one "achieve" function
# should be used in a single pipeline; otherwise, the second one would
# overwrite the objective of the first one.


# Set up a MIPContext with based on an additive valuation profile. Initializes a
# JuMP model with constraints ensuring a valid allocation. These constraints can
# be modified by suppling lower and upper limits to the number of items in a
# bundle, or the number of owners an item may have. These arguments can be
# `nothing` (no limit), an integer limit, or a vector of integers, with one
# limit per agent or item.
function init_mip(V::Additive, solver;
        min_bundle=nothing,   # Agents may receive nothing
        max_bundle=nothing,   # Agents may receive all items
        min_owners=1,         # Items have at least one owner
        max_owners=1)         # Items have at most one owner

    model = Model(solver)

    N = agents(V)
    M = items(V)

    @variable(model, A[N, M], binary=true)

    isnothing(min_bundle) || broadcast(N, min_bundle) do i, lo
        @constraint(model, sum(A[i, g] for g in M) >= lo)
    end

    isnothing(max_bundle) || broadcast(N, max_bundle) do i, hi
        @constraint(model, sum(A[i, g] for g in M) <= hi)
    end

    isnothing(min_owners) || broadcast(M, min_owners) do g, lo
        @constraint(model, sum(A[i, g] for i in N) >= lo)
    end

    isnothing(max_owners) || broadcast(M, max_owners) do g, hi
        @constraint(model, sum(A[i, g] for i in N) <= hi)
    end

    return MIPContext(V, A, model, solver)

end


# See `mms` and `alloc_mms` for special cases.
const MIP_LIMIT_DOC = """
Lower and upper limits on the size of each bundle and the number of owners for
each item may be supplied using the keyword arguments `min_bundle`,
`max_bundle`, `min_owners` and `max_owners`, the latter two of which default to
`1`. If one of these is `nothing`, the limit is simply absent. Otherwise, the
argument is broadcast to the appropriate size.
"""


# Shortcut, so you can use a matrix instead of an Additive valuation profile,
# including in the alloc_... functions that rely on init_mip.
init_mip(V::Matrix, solver; kwds...) = init_mip(Additive(V), solver; kwds...)


# Solves the MIP model and constructs the actual Allocation object in the
# MIPContext.
# ϵ: https://www.gurobi.com/documentation/9.1/refman/intfeastol.html
function solve_mip(ctx; ϵ=1e-5, check=nothing)

    if !isempty(ctx.callbacks)
        set_attribute(ctx.model, MOI.LazyConstraintCallback(),
            cb_data -> lazy_constraint_callback(ctx, cb_data))
    end

    if isempty(ctx.objectives)
        optimize!(ctx.model)
    else
        # Ignoring the returned constraints -- not doing any cleanup
        lex_optimize!(ctx.model, ctx.objectives)
    end

    st = termination_status(ctx.model)
    @assert conf.MIP_SUCCESS === nothing || 
            st in conf.MIP_SUCCESS "solver termination status is $st"

    V = ctx.profile
    ctx.alloc = Allocation(na(V), ni(V))

    for i in agents(V), g in items(V)
        val = JuMP.value(ctx.alloc_var[i, g])
        @assert val ≤ ϵ || val ≥ 1 - ϵ
        val ≥ 1.0 - ϵ && give!(ctx.alloc, i, g)
    end

    isnothing(check) || check(ctx.alloc)

    return ctx

end


# The "master" lazy constraint callback method that runs each registered
# callback.
function lazy_constraint_callback(ctx, cb_data)
    if callback_node_status(cb_data, ctx.model) == MOI.CALLBACK_NODE_STATUS_INTEGER
        for f in ctx.callbacks
            f(ctx, cb_data)
        end
    end
end


## Objective-modifying pipeline steps (achieve_...)


# Set up objective and constraints to make sure the JuMP model produces an MNW
# allocation, using a slightly adapted version of the approach of Caragiannis
# et al. (https://doi.org/10.1145/3355902).
achieve_mnw(mnw_warn) = function(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    @assert isintegral(V)
    @assert isnonnegative(V)

    N = [i for (i, g) in bipartite_matching(matrix(V))]

    M = items(V)

    v_max = Float64(maximum(value(V, i, M) for i in N))

    mnw_prec = true

    for (k, name) in [1 => "PO", 2 => "EF1", length(N) => "MNW"]
        if log(v_max^k) - log(v_max^k - 1) == 0.0
            if name == "MNW"
                mnw_prec = false
                !mnw_warn && continue
            end
            @warn("Precision possibly insufficient to guarantee $name")
        end
    end

    ctx.res = (ctx.res..., mnw_prec = mnw_prec)

    @variable(model, W[N])

    @objective(model, Max, sum(W))

    A = ctx.alloc_var

    for i in N, k = 1:2:v_max
        @constraint(model, W[i] <=
                log(k) + (log(k + 1) - log(k)) * (value(V, i, A) - k))
    end

    for i in N
        @constraint(model, value(V, i, A) >= 1)
    end

    return ctx

end


# Set up objective and constraints to make sure the JuMP model produces an
# egalitarian/maximin allocation.
achieve_mm(cutoff=nothing, ignored_agents=[]) = function(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    @variable(model, v_min)

    for i in setdiff(N, ignored_agents)
        @constraint(model, v_min <= value(V, i, A))
    end

    isnothing(cutoff) || @constraint(model, v_min <= cutoff)

    if !issubset(N, ignored_agents)
        @objective(model, Max, v_min)
    end

    return ctx

end


# Set up objective and constraints to make sure the JuMP model produces an
# lexicographic maximin/leximin allocation.
function achieve_lmm(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model
    N = agents(V)

    # Negated, because we're using a leximax formulation:
    funcs = [-value(V, i, A) for i in N]
    objectives = leximax_os06_1!(model, funcs)

    ctx.objectives = [(MOI.MIN_SENSE, obj) for obj in objectives]

    return ctx

end


# Set up objective and constraints to make sure the JuMP model produces a
# minimum ordered weighted average, by utility rank, using wt as the weight
# function (cf., Lesca & Perny, 2010)
achieve_ggi(wt) = function(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    N, n, M = agents(V), na(V), items(V)

    @variable(model, r[N])
    @variable(model, b[N, N] >= 0)

    omega = [i == n ? wt(i, n) : wt(i, n) - wt(i + 1, n) for i in N]

    @objective(model, Max,
        sum(omega[j] * (j * r[j] - sum(b[i, j] for i in N)) for j in N))

    for i in N, j in N
        @constraint(model,
            r[j] - b[i, j] <= value(V, i, A))
    end

    return ctx

end


# Set up objective for maximum utilitarian welfare (MUW).
function achieve_muw(ctx)
    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    @objective(model, Max, utility(V, A))

    return ctx
end


## Objective-preserving pipeline steps (enforce_...)

# How `enforce` works depends on whether a constraint is symmetric or not.
# Symmetric constraints can simply implement `enforce(C)`, and be done with it.
# Asymmetric constraints, however, should implement `enforce(C, i, j)` instead,
# which enforces the individual bundle constraint of agent `i` on the bundle of
# agent `j`. This is used in `alloc_mms`, where the normal behavior (each agent
# has her own bundle constraint) is used on the allocation itself, but where
# each agent's constraint us duplicated along with her valuation function when
# finding here MMS (with `mms`).


# Enforce no constraints on the JuMP model.
enforce(C::Nothing) = identity


# Default: Each agent has her own bundle constraint. Symmetric constraints can
# simply override this implementation directly.
enforce(C::Constraint) = function(ctx)
    reduce(|>, (enforce(C, i, i) for i in agents(ctx)), init=ctx)
end


# Need not be implemented directly by any constraint.
enforce(C, i) = _enforce(C, i, Symmetry(C))
_enforce(C, i, ::Symmetric) = enforce(C)
_enforce(C, i, ::Asymmetric) = function(ctx)
    reduce(|>, (enforce(C, i, j) for j in agents(ctx)), init=ctx)
end


# Enforce a tuple of constraints.
enforce(C::Constraints, args...) = function(ctx)
    reduce(|>, (enforce(Cᵢ, args...) for Cᵢ in C.parts), init=ctx)
end


# Enforce cardinality constraints on the JuMP model.
enforce(C::Counts) = function(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    for S in C, i in agents(V)
        @constraint(model, sum(A[i, g] for g in S) <= S.threshold)
    end

    return ctx

end


# Enforce item conflict constraints on the JuMP model.
enforce(C::Conflicts) = function(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    G = C.graph

    @assert nv(G) == ni(V)

    for e in edges(G)
        g, h = src(e), dst(e)
        for i in agents(V)
            @constraint(model, A[i, g] + A[i, h] <= 1)
        end
    end

    return ctx

end


# Internal. Used by `enforce` for `Required`, `Forbidden` and `Permitted`.
function fix_match_vars(ctx, C, i, j, val=1, which=bundle)

    A₀, A, model = C.alloc, ctx.alloc_var, ctx.model

    # The items affected by agent i's part of the constraint
    affected = which(A₀, i)

    # Apply the constraint to the bundle of agent j
    for g in affected
        fix(A[j, g], val)
    end

    return ctx

end


# Enforce inclusion constraints on the JuMP model.
enforce(C::Required, i, j) = function(ctx)
    fix_match_vars(ctx, C, i, j)
end


# Enforce exclusion constraints on the JuMP model.
enforce(C::Forbidden, i, j) = function(ctx)
    fix_match_vars(ctx, C, i, j, 0)
end


# Enforce permission constraints (i.e., exclusion constraints with the
# complement) on the JuMP model.
enforce(C::Permitted, i, j) = function(ctx)
    which(A, i) = setdiff(items(A), bundle(A, i))
    fix_match_vars(ctx, C, i, j, 0, which)
end


# Add the constraint that any bundle must contain at most r items, where r
# is the rank of the matroid, since any independent set of a matroid will
# have at most r items.
function matroid_initial_constraint(ctx, M, i)
    A = ctx.alloc_var
    E = ground_set(M)
    r = rank(M)
    @constraint(ctx.model, sum(A[i, g] for g in E) <= r)
end


# Check if agent `i` has a dependent bundle in the matroid `M`, and submit a
# constraint to the model if necessary.
function matroid_fix_constraint(ctx, cb_data, M, i)
    V, A, model = ctx.profile, ctx.alloc_var, ctx.model
    ϵ = 1e-5
    bundle = Set()

    for g in items(V)
        val = callback_value(cb_data, A[i, g])
        @assert val <= ϵ || val >= 1 - ϵ
        val >= 1.0 - ϵ && push!(bundle, g)
    end

    if !is_indep(M, bundle)
        bundle_rank = rank(M, bundle)
        con = @build_constraint(sum(A[i, g] for g in bundle) <= bundle_rank)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
end


# Enforce a symmetric matroid constraint on the JuMP model (using a lazy
# callback).
enforce(C::MatroidConstraint) = function (ctx)
    callback = function (ctx, cb_data)
        V = ctx.profile
        M = C.matroid
        for i in agents(V)
            matroid_fix_constraint(ctx, cb_data, M, i)
        end
    end

    push!(ctx.callbacks, callback)

    V = ctx.profile
    M = C.matroid
    for i in agents(V)
        matroid_initial_constraint(ctx, M, i)
    end

    return ctx
end


# Enforce an asymmetric matroid constraint on the JuMP model (using lazy
# callbacks).
enforce(C::MatroidConstraints, i, j) = function (ctx)
    callback = function (ctx, cb_data)
        M = C.matroids[i]
        matroid_fix_constraint(ctx, cb_data, M, j)
    end

    push!(ctx.callbacks, callback)

    M = C.matroids[i]
    matroid_initial_constraint(ctx, M, j)

    return ctx
end


# Enforce envy-freeness up to a single object (EF1) on the JuMP model.
function enforce_ef1(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

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
        vii = value(V, i, A)

        # Agent i's value for j's bundle, without dropped item
        vij1 = sum(value(V, i, g) * (A[j, g] - D[i, j, g]) for g in M)

        # No envy, once an item is (possibly) dropped:
        @constraint(model, vii >= vij1)

    end

    return ctx

end


# Enforce envy-freeness up to any object (EFX) on the JuMP model.
function enforce_efx(ctx)

    V, A, model = ctx.profile, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    for i in N, j in N

        i == j && continue

        # Agent i's value for her own bundle
        vii = value(V, i, A)

        # Value large enough to skip the constraint:
        huge = value(V, i, M)

        # For each possibly dropped item d ...
        for d in M

            # Agent i's value for j's bundle, without d
            vijx = sum(value(V, i, g) * A[j, g] for g in M if g ≠ d)

            # If d isn't allocated to j, we skip this constraint:
            skip = huge * (1 - A[j, d])

            # No envy, once d is dropped (or if d isn't assigned to j):
            @constraint(model, vii + skip >= vijx)

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
    alloc_ef1(V, C; solver=conf.MIP_SOLVER, kwds...)

Create an `Allocation` that is envy-free up to one item (EF1), based on the
valuation profile `V`, possibly subject to the constraints given by the
`Constraint` object `C`. The solution is found using a straightforward
mixed-integer program, and is most suitable for constraints where no specialized
algorithm exists. For example, without constraints, a straightforward round
robin picking sequence yields EF1, and a similar strategy works for cardinality
constraints. (It is still possible to use this function without constraints, by
explicitly supplying `nothing` for the constraint argument `C`.) The return
value is a named tuple with the fields `alloc` (the `Allocation`) and `model`
(the JuMP model used in the computation).

$MIP_LIMIT_DOC

Note that for some constraints, there may not *be* an EF1 allocation, in which
case the function will fail with an exception.
"""
function alloc_ef1(V, C; solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    enforce_ef1 |>
    enforce(C) |>
    solve_mip |>
    alloc_result

end


"""
    alloc_efx(V[, C]; solver=conf.MIP_SOLVER, kwds...)

Create an `Allocation` that is envy-free up to any item (EFX), based on the
valuation profile `V`, possibly subject to the constraints given by the
`Constraint` object `C`. The solution is found using a straightforward
mixed-integer program. The return value is a named tuple with the fields `alloc`
(the `Allocation`) and `model` (the JuMP model used in the computation).

$MIP_LIMIT_DOC

Note that while some constraints may prevent an exact EFX allocation, it is
currently (Mar 2021) an open question whether EFX always exists in the
unconstrained case (see, e.g., [Improving EFX Guarantees through Rainbow Cycle
Number](https://arxiv.org/abs/2103.01628) by Chaudhury et al.).
"""
function alloc_efx(V, C=nothing; solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    enforce_efx |>
    enforce(C) |>
    solve_mip |>
    alloc_result

end


# Extract the allocation and the MNW value (excluding agents with a utility of
# zero) at the end of the pipeline. Strictly speaking, we needn't include the
# mnw field here, as it's a separate calculations; however, it's convenient if
# you supply a matrix as the argument for alloc_mnw function, since we have
# access to the valuation profile here. (Also, the cost of nash_welfare will
# generally be very low, compared to the actual MIP solving.)
function mnw_result(ctx)
    V, A = ctx.profile, ctx.alloc
    return (alloc=A, model=ctx.model, mnw=nash_welfare(V, A), ctx.res...)
end


"""
    alloc_mnw(V[, C]; mnw_warn=false, solver=conf.MIP_SOLVER, kwds...)

Create an `Allocation` attaining maximum Nash welfare (MNW), based on the
valuation profile `V`, possibly subject to the constraints given by the
`Constraint` object `C`. The solution is found using the approach of Caragiannis
et al. in their 2019 paper [The Unreasonable Fairness of Maximum Nash
Welfare](https://doi.org/10.1145/3355902), with two minor modifications:

1. Rather than hard-coding a maximum valuation (arising from the assumption that
   the values of each agent sum to 1000), this maximum is extracted from `V`;
   and

2. Extra constraints are permitted (through the object `C`), possibly lowering
   the attainable MNW.

Because of how the integer program is constructed, it may be affected by
precision effects, where a high number of agents can make it impossible to
guarantee Pareto optimalty (PO), EF1 or MNW. If the precision is too low, the
appropriate warning will be issued, but the computation is not halted. Note that
these warnings are quite conservative (see note below). This is particularly
true of the one for MNW, which is disabled by default, in part because of its
sensitivity, and in part because it will generally be useful to find solutions
that satisfy PO and EF1, even if it may not be exactly MNW. The MNW warning can
be enabled by setting the `mnw_warn` keyword to `true`.

!!! note

    The warnings are based on the lower bounds described by Caragiannis et al.
    On the one hand, the bound is only used to test whether current
    floating-point precision is sufficient; any tolerance or gap used by the
    solver is not used, which might in principle mean that false negatives are
    possible. On the other hand, these bounds, especially the one for exact MNW,
    may in practice be quite loose, with small variations in agent utilities
    leading to large changes in objective value, unless the changes are finely
    tuned to cancel out.

The return value is a named tuple with fields `alloc` (the `Allocation`),
`mnw` (the achieved Nash welfare for the agents with nonzero utility),
`mnw_prec` (whether or not there was enough precision to satisfy the lower bound
guaranteeing exact MNW) and `model` (the JuMP model used in the computation).

$MIP_LIMIT_DOC
"""
function alloc_mnw(V, C=nothing; mnw_warn=false, solver=conf.MIP_SOLVER,
        kwds...)

    init_mip(V, solver; kwds...) |>
    achieve_mnw(mnw_warn) |>
    enforce(C) |>
    solve_mip |>
    mnw_result

end


"""
    alloc_mnw_ef1(V, C; mnw_warn=true, solver=conf.MIP_SOLVER, kwds...)

Equivalent to `alloc_mnw`, except that EF1 is enforced. Without any added
constraints, MNW implies EF1, so this function is not needed in that case.
Therefore the argument `C` is not optional.

$MIP_LIMIT_DOC
"""
function alloc_mnw_ef1(V, C; mnw_warn=true, solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    achieve_mnw(mnw_warn) |>
    enforce_ef1 |>
    enforce(C) |>
    solve_mip |>
    mnw_result

end


# Extract the allocation, model and maximin value at the end of the pipeline.
function mm_result(ctx)

    # If there is no objective, all agents were ignored, and the maximin value
    # is unbounded:
    if objective_sense(ctx.model) == MOI.FEASIBILITY_SENSE
        mm = Inf
    else
        mm = objective_value(ctx.model)
    end

    return (alloc=ctx.alloc, model=ctx.model, mm=mm, ctx.res...)

end


"""
    alloc_mm(V[, C]; cutoff=nothing, ignored_agents=[],
        solver=conf.MIP_SOLVER, kwds...)

Create an egalitarian or maximin `Allocation`, i.e., one where the minimum
bundle value is maximized. The `cutoff`, if any, is a level at which we are
satisfied, i.e., any allocation where all agents attain this value is
acceptable. The return value is a named tuple with fields `alloc` (the
`Allocation`), `mm` (the lowest agent utility) and `model` (the JuMP model
used in the computation).

The `ignored_agents` argument indicates agents that should be ignored when
maximizing the minimum. These agents may still receive items, and will
participate in forming a feasible allocation (possibly with respect to some
constraint `C`); they are only ignored in the objective.

!!! note

    Most users will probably not need `ignored_agents`. Its primary use is as
    part of `alloc_mms`, for ignoring agents with an MMS of zero, whose α is
    unbounded.

$MIP_LIMIT_DOC
"""
function alloc_mm(V, C=nothing; cutoff=nothing, ignored_agents=[],
        solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    achieve_mm(cutoff, ignored_agents) |>
    enforce(C) |>
    solve_mip |>
    mm_result

end


"""
    alloc_lmm(V[, C]; solver=conf.MIP_SOLVER, kwds...)

Create a lexicographic maximin (leximin) `Allocation`, i.e., one where the
lowest bundle value is maximized, and subject to that, the second lowest is
maximized, etc. The return value is a named tuple with fields `alloc` (the
`Allocation`) and `model` (the JuMP model used in the computation).

The method used for leximin optimization is that of Ogryczak and Śliwiński (["On
Direct Methods for Lexicographic Min-Max
Optimization"](https://doi.org/10.1007/11751595_85), 2006).

$MIP_LIMIT_DOC
"""
function alloc_lmm(V, C=nothing; solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    achieve_lmm |>
    enforce(C) |>
    solve_mip |>
    alloc_result

end


# Acts as if agent `i`'s part of constraint `C` applies to all bundles. If `C`
# is `Symmetric`, this has no effect.
struct SymmetrizedConstraint{T <: Union{Constraint,Nothing}} <: Constraint
    C::T
    i::Int
end
enforce(C::SymmetrizedConstraint) = enforce(C.C, C.i)


# Get the limit (from `min_bundle` or `max_bundle`) of a given agent
function get_limit(lim, i)
    (isnothing(lim) || length(lim) == 1) && return lim
    return lim[i]
end


"""
    mms(V::Additive, i[, C]; solver=conf.MIP_SOLVER, kwds...)

Determine the maximin share of agent `i`, i.e., the bundle value she is
guaranteed to attain if she partitions the items and the other agents choose
their bundles. Useful, e.g., as a point of reference when determining the
empirical approximation ratios of approximate MMS allocation algorithms. Also
used as a subroutine in `alloc_mms`. The return value is a named tuple with the
fields `mms` (the maximin share of agent `i`) and `model` (the JuMP model used
in the computation).

$MIP_LIMIT_DOC

If a constraint `C` is supplied, and this is asymmetric (i.e., different for the
different agents), agent `i`'s version is enforced on every bundle to determine
which partitions are feasible.
"""
function mms(V::Additive, i, C=nothing; solver=conf.MIP_SOLVER, kwds...)

    # Let all agents be clones of agent i
    Vᵢ = Additive([value(V, i, g) for _ in agents(V), g in items(V)])
    Cᵢ = SymmetrizedConstraint(C, i)
    minbᵢ = get_limit(get(kwds, :min_bundle, nothing), i)
    maxbᵢ = get_limit(get(kwds, :max_bundle, nothing), i)

    # maximin in this scenario is MMS for agent i
    res = alloc_mm(Vᵢ, Cᵢ; solver=solver, kwds...,
                   min_bundle=minbᵢ, max_bundle=maxbᵢ)

    return (mms=res.mm, model=res.model)

end


"""
    mms_alpha(V, A, mmss)

Utility function to find the fraction of the maximin share guarantee attained by
the allocation `A`, under the valuation profile `V`, where `mmss[i]` is the MMS
of agent `i`. This makes it possible, for example, to use the `mmss` field from
the result of `alloc_mms` to find the MMS approximation provided by an
allocation constructed by other means. For example:

    mmss = alloc_mms(V).mmss
    A = alloc_rand(V).alloc
    alpha = mms_alpha(V, A, mmss)

If all agents have an MMS of zero, `alpha` will be unbounded, represented by the
value `Inf`. This is the case even if one or more agents get a value of `0`.
"""
function mms_alpha(V, A, mmss)
    vals = [value(V, i, A) for i in agents(A)]
    frac = vals ./ mmss
    frac[isnan.(frac)] .= Inf # we want x/0 to yield Inf even when x == 0
    return minimum(frac)
end


"""
    alloc_mms(V[, C]; cutoff=false, solver=conf.MIP_SOLVER,
              mms_kwds=(), kwds...)

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

If all agents have an MMS of zero, `alpha` will be unbounded, represented by the
value `Inf`. This is the case even if one or more agents get a value of `0`.

$MIP_LIMIT_DOC

When finding an MMS partition for each individual agent, the constraint and
options act a little differently than when finding the actual allocation. If an
asymmetric `Constraint` `C` is supplied (i.e., one that is different for the
different agents), agent `i`'s version is enforced on every bundle to determine
which partitions are feasible.

The same holds for `min_bundle` and `max_bundle`. Apart from those two, the
keywords `kwds` are also used for the MMS partitions.

!!! tip

    In some cases, using asymmetric constraints might lead to a situation where
    a feasible allocation exists, but for some agent, an MMS partition does not.
    Then the agent's maximin share is undefined! One way around this is to relax
    the definition a bit, and to permit charity when finding the MMS partitions.
    The agent will still try to maximize the mininum value of any bundle in this
    partition, but is not required to allocate every object. This can be
    achieved as follows:

    ```julia
    alloc_mms(V, C, mms_kwds=(min_owner=0,))
    ```

    You should verify that the strategy makes sense for your application,
    however! In some cases, it might not be necessary, and would just needlessly
    inflate maximin shares. Also, it is only guaranteed to work for constraints
    where an empty bundle is feasible; it might fail for [`Required`](@ref), for
    example. Indeed, for some constraints, one could argue that it might not be
    advisable to use the MMS criterion to begin with.

"""
function alloc_mms(V::Additive, C=nothing; cutoff=false, solver=conf.MIP_SOLVER,
        mms_kwds=(), kwds...)

    N, M = agents(V), items(V)

    X = zeros(na(V), ni(V))

    ress = [mms(V, i, C; solver=solver, kwds..., mms_kwds...) for i in N]

    # Individual maximin shares -- also included in the result
    mmss = [res.mms for res in ress]

    mms_models = [res.model for res in ress]

    for (i, μ) in enumerate(mmss), g in M
        X[i, g] = value(V, i, g) / μ
    end

    max_alpha = cutoff ? 1.0 : nothing

    res = alloc_mm(Additive(X), C;
        cutoff = max_alpha,
        ignored_agents = N[iszero.(mmss)],
        solver = solver,
        kwds...)

    return (alloc       = res.alloc,
            model       = res.model,
            mms_models  = mms_models,
            alpha       = res.mm,
            mmss        = mmss)

end


alloc_mms(V::Matrix, C=nothing; kwds...) = alloc_mms(Additive(V), C; kwds...)


# Extract the allocation at the end of the pipeline.
function ggi_result(ctx)
    return (alloc=ctx.alloc, model=ctx.model, ctx.res...)
end


"""
    wt_gini(i, n)

The (unnormalized) weights used in the ordered weighted average in the Gini
social-evaluation function, where the utility of the `i`th agent, ordered by
increasing utility, is given weight ``2(n - i) + 1``. (The normalized weights
yielding the original Gini social-evaluation function are divided by ``n^2``,
but this makes no difference to the optimization problem.)
"""
wt_gini(i, n) = 2(n - i) + 1


"""
    alloc_ggi(V[, C]; wt=wt_gini, solver=conf.MIP_SOLVER, kwds...)

Maximizes a generalized Gini index (GGI), also known as a generalized Gini
social-evaluation functions. The function being maximized is an ordered weighted
average (OWA) of agent utilities, utilities, where the weight is based on
utility rank `i`, from the least happy (`1`) to the most happy (`n`),
parameterized by the function `wt(i, n)`. It is generally assumed that the
weights are nondecreasing in `i`. Note that there is no need to use normalized
weights (i.e., to produce a weighted average, despite the term OWA), as is often
the case when such measures are used to measure *in*equality (e.g., by
subtracting the OWA from an ordinary average, cf. [Generalized gini inequality
indices](https://doi.org/10.1016/0165-4896(81)90018-4) by John A. Weymark).

The default `wt_gini` gives the (non-normalized) weights of the original Gini
social-evaluation. Two other notable cases for `wt` are `(i, _) -> i == 1`,
which yields a maximin allocation, and `(i, _) -> 1`, which yields a purely
utilitarian allocation (with no consideration for fairness). The solution method
used is based on that of Lesca and Perny (linear formulation ``\\Pi'_W``) in
their paper 2010 paper [“LP Solvable Models for Multiagent Fair Allocation
Problems”](https://doi.org/10.3233/978-1-60750-606-5-393). The return value is a
named tuple with the fields `alloc` (the `Allocation` that has been produced)
and `model` (the JuMP model used in the computation).

$MIP_LIMIT_DOC

In the original inequality measures, the mean agent utility is included as a
normalizing term, which is harmless for the case of identical valuations
functions (and when looking at, say, the distribution of incomes), but when
valuations differ, this mean will vary with the allocations. As pointed out by
Lesca and Perny, such a measure is not monotone with Pareto dominance -- the
optimization will tend to drive the mean utility *down*. Therefore only the term
measuring (in)equality (i.e., the ordered weighted sum of agent utilities) is
used.
"""
function alloc_ggi(V, C=nothing; wt=wt_gini, solver=conf.MIP_SOLVER, kwds...)

    init_mip(V, solver; kwds...) |>
    achieve_ggi(wt) |>
    enforce(C) |>
    solve_mip |>
    ggi_result

end


"""
    rand_priority_profile(n, m; rng=default_rng())

Generate a random valuation profile, based on a random prioritization of the
items. To ensure that the item priorities will overrule any other
considerations, the valuation profile is filled using powers of 2. Therefore,
this function does not work for values of `m` ≥ 63.
"""
function rand_priority_profile(n, m; rng=default_rng())

    m < 63 || throw(DomainError(m, "number of items ≥ 63 is not supported"))

    # Random permutation of the items for priority
    π = randperm(rng, m)
    # Fill the valuation matrix with ones, so that no item is unwanted
    X = ones(n, m)

    # Encode the item priorities in an additive profile matrix
    for g in 1:m
        i = rand(rng, 1:n)
        X[i, g] = 2 ^ π[g]
    end

    return Profile(X)

end


function rand_mip_result(ctx)
    return (profile=ctx.profile, alloc=ctx.alloc, model=ctx.model, ctx.res...)
end


"""
    alloc_rand_mip(V[, C]; solver=conf.MIP_SOLVER, rng=default_rng(),
                   kwds...)

Allocate items to agents randomly, in a similar manner to [`alloc_rand`](@ref).
The implementation is inspired by the [randomized coloring procedure with
symmetry-breaking](https://doi.org/10.1007/978-3-540-70575-8_26) of Pemmaraju
and Srinivasan, and is generalized to work with any `Constraint`.

The profile `V` is not used directly, other than to determine the number of
agents and items.

A random priority profile is generated using [`rand_priority_profile`](@ref),
such that any constraints will remove the items with lowest priority first.
To ensure that the MIP respects this priority profile, the utilitarian welfare
of the allocation is maximized.

$MIP_LIMIT_DOC
"""
function alloc_rand_mip(V, C=nothing; solver=conf.MIP_SOLVER, rng=default_rng(),
        kwds...)

    init_mip(rand_priority_profile(na(V), ni(V), rng=rng), solver; kwds...) |>
    achieve_muw |>
    enforce(C) |>
    solve_mip |>
    rand_mip_result

end
