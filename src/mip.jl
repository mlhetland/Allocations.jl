# Algorithms that share some allocation-specific infrastructure based on
# mixed-integer programming using JuMP.


# A shared context for the MIP pipeline. The JuMP variable representing the
# allocation (as an n-by-m matrix) is kept in alloc_var.
mutable struct MIPContext{V <: Additive, M, A, S}
    valuation::V
    alloc_var::A
    model::M
    solver::S
    alloc::Union{Allocation, Nothing}
end
MIPContext(v, a, m, s) = MIPContext(v, a, m, s, nothing)


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

    @assert termination_status(ctx.model) == MOI.OPTIMAL

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
function achieve_mnw(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    @assert isintegral(V)
    @assert isnonnegative(V)

    positive = [value(V, i, j) > 0 for i in agents(V), j in items(V)]

    matching = bipartite_matching(positive, ctx.solver)

    N = findall(vec(any(matching, dims=2)))

    M = items(V)

    v_max = maximum(value(V, i, M) for i in N)

    for (k, name) in [1 => "PO", 2 => "EF1", length(N) => "MNW"]
        if log(v_max^k) - log(v_max^k - 1) == 0.0
            @warn("Precision insufficient to guarantee $name")
        end
    end

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
function achieve_mm(ctx)

    V, A, model = ctx.valuation, ctx.alloc_var, ctx.model

    N, M = agents(V), items(V)

    @variable(model, v_min)

    for i in N
        @constraint(model, v_min <= sum(A[i, g] * value(V, i, g) for g in M))
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


# Actual allocation methods


# Extract the allocation and the MNW value (excluding agents with a utility of
# zero) at the end of the pipeline. Strictly speaking, we needn't include the
# mnw field here, as it's a separate calculations; however, it's convenient if
# you supply a matrix as the argument for alloc_mnw function, since we have
# access to the valuation here. (Also, the cost of nash_welfare will generally
# be very low, compared to the actual MIP solving.)
function mnw_result(ctx)
    V, A = ctx.valuation, ctx.alloc
    return (alloc=A, mnw=nash_welfare(V, A))
end


"""
    alloc_mnw(V[, C]; solver=conf.MIP_SOLVER)

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
guarantee Pareto optimalty, EF1 or MNW respectively). If the precision is too
low, the appropriate warning will be issued, but the computation is not
halted.

The return value is a named tuple with fields `alloc` (the `Allocation`) and
`mnw` (the achieved Nash welfare for the agents with nonzero utility).
"""
function alloc_mnw(V, C=nothing; solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mnw |>
    enforce(C) |>
    solve_mip |>
    mnw_result

end


# Extract the allocation and the maximin value at the end of the pipeline.
function mm_result(ctx)
    return (alloc=ctx.alloc, mm=objective_value(ctx.model))
end


"""
    alloc_mm(V[, C]; solver=conf.MIP_SOLVER)

Create an egalitarion or maximin `Allocation`, i.e., one where the minimum
bundle value is maximized. The return value is a named tuple with fields
`alloc` (the `Allocation`) and `mm` (the lowest agent utility).
"""
function alloc_mm(V, C=nothing; solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mm |>
    enforce(C) |>
    solve_mip |>
    mm_result

end


"""
    alloc_mms(V[, C]; solver=conf.MIP_SOLVER)

Find an MMS allocation, i.e., one that satisfies the *maximin share
guarantee*, where each agent gets a bundle it weakly prefers to its maximin
share (introduced by Budish, in his 2011 paper [The Combinatorial Assignment
Problem: Approximate Competitive Equilibrium from Equal
Incomes](https://doi.org/10.1086/664613)). The return value is a named tuple
with fields `alloc` (the `Allocation`) and `alpha`, the lowest fraction of MMS
that any agent achieves (is at least 1 exactly when the allocation is MMS).
"""
function alloc_mms(V::Additive, C=nothing; solver=conf.MIP_SOLVER)

    N, M = agents(V), items(V)

    X = zeros(na(V), ni(V))

    for i in N

        # Let all agents be clones of agent i
        Vi = Additive([value(V, i, g) for _ in N, g in M])

        # maximin in this scenario is MMS for agent i
        res = alloc_mm(Vi, C, solver=solver)

        # Scale agent i's values by agent i's MMS
        for g in M
            X[i, g] = value(V, i, g) / res.mm
        end

    end

    # maximin with scaled values is as close to the MMS guarantee as possible
    res = alloc_mm(Additive(X), C, solver=solver)

    return (alloc=res.alloc, alpha=res.mm)

end


alloc_mms(V::Matrix, C=nothing; solver=conf.MIP_SOLVER) =
    alloc_mms(Additive(V), C, solver=solver)


# Extract the allocation at the end of the pipeline.
function mgg_result(ctx)
    return (alloc=ctx.alloc,)
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
The return value is a named tuple with the field `alloc`, the `Allocation`
that has been produced.
"""
function alloc_mgg(V, C=nothing; wt=wt_gini, solver=conf.MIP_SOLVER)

    init_mip(V, solver) |>
    achieve_mgg(wt) |>
    enforce(C) |>
    solve_mip |>
    mgg_result

end
