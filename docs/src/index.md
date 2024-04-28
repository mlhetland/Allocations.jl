# Allocations.jl

```@meta
CurrentModule = Allocations
```

The Allocations package deals with the fair allocation of indivisible items
to a set of agents. An instance of this problem consists of:

- A set $N$ of $n$ *agents* and a set $M$ of $m$ *items*;
- One *valuation function* $v_i:2^M\to\R$ for each agent $i\in N$;

For simplicity, we may simply let $N=\{1,\dots,n\}$ and $M=\{1,\dots,m\}$.[^1]

[^1]: The latter is less common, presumably because the two sets then intersect. An alternative is to let $M$ be a set of opaque objects $g_j$, for $j=1\dots m$.

The goal is to find an *allocation* $A=(A_1,\dots,A_n)$, which gives a *bundle*
$A_i\subseteq M$ to each agent, and which satisfies certain fairness criteria.
The allocation is usually required to form a partition of $M$, though in some
scenarios, one may deviate from this.

In `Allocations.jl`, an instance is represented by the *valuation profile*
$V=\{v_i:i\in N\}$, in the form of a `Profile` object. Given a profile `V`, the
various parts of the instance may be accessed as follows:

- `items(V)`: The item set $N$, as an iterable;
- `agents(V)`: The agent set $M$, as an iterable;
- `value(V, i, S)`: The value agent $i$ assigns to the set $S\subseteq
  M$, i.e., $v_i(S)$.

`S` may be any iterable. If it a single item (i.e., `Int`) `g`, it is
interpreted as a singleton, `[g]` (though usually handled more efficiently).
The number of agents and items, respectively, are found using `na(V)` and
`ni(V)`

In addition to this basic setting, instances may be constrained, by supplying
some *constraint object*, describing which bundles are feasible. For example, if
one is looking for bundles that are *connected* in some sense, the constraint
object will typically be a graph that defines the connectivity relation, etc.

Allocations are represented by `Allocation` objects. Given an allocation `A`,
the bundle of agent `i` is found using `bundle(A, i)`.

!!! tip

    For more on this topic, see, e.g., the [Wikipedia entry on fair item
    allocation](https://en.wikipedia.org/wiki/Fair_item_allocation), or the
    surveys by [Amanatidis et al.](https://arxiv.org/abs/2208.08782) and
    [Suksompong](https://doi.org/10.1145/3505156.3505162), on the unconstrained
    and constrained versions of the problem, respectively.

## Installation

To install the package, you can simply import it in the [Julia
REPL](https://docs.julialang.org/en/v1/stdlib/REPL/):

```jldoctest intro
julia> using Allocations
```

Press enter at the resulting prompt to install both the package and its
dependencies.

To install a more recent version than the released one, you can use the package
manager directly. In the Julia REPL, press `]` to enter the `Pkg` REPL, and then
add the package directly from the source:

```julia-repl
pkg> add https://github.com/mlhetland/Allocations.jl
```

You can then import the module as before.

## Basic use

To specify an allocation problem instance, create a valuation profile:

```jldoctest intro
julia> V = Profile([1 2 3; 2 3 1])
Additive{Matrix{Int64}} with 2 agents and 3 items:
 1  2  3
 2  3  1
```

`Profile` is an abstract class, and `Profile(X::Matrix)` is an alias for
`Additive(X)`. Once you have a valuation profile, you can use an allocation
function (ones called `alloc_...`), e.g., for finding a maximum Nash welfare
(MNW) allocation:

```jldoctest intro
julia> res = alloc_mnw(V);
```

Note that the first time you call an allocation function, it may take some time
to finish, because there's quite a lot of compilation going on behind the
scenes. From then on, in the same REPL session, there will be much less
overhead.

These functions take a `Profile` as input and return a named tuple with the
field `alloc` referring to an `Allocation`:

```jldoctest intro
julia> A = res.alloc
Allocation with 2 agents and 3 items:
  1 => {3}
  2 => {1, 2}
```

The bundle of each agent is available through the `bundle` function:

```jldoctest intro
julia> bundle(A, 2)
Set{Int64} with 2 elements:
  2
  1
```

Bundles should not be modified directly, as the `Allocation` also maintains an
inverse mapping, from items to agents. Rather, use the `give!` and `deny!`
functions.

Some allocation functions may produce other results as well, such as properties
of the allocation that are naturally computed as part of the allocation process.
For the MNW case, the objective value (the Nash welfare, which is being
maximized) is available as `mnw`:

```jldoctest intro
julia> res.mnw
15.0
```

The allocation functions also permit a matrix argument as a shortcut, implicitly
creating an `Additive`. For example, you can find a maximin share (MMS)
allocation as follows:

```jldoctest intro
julia> alloc_mms([1 1 2 3; 2 1 2 3]).alloc
Allocation with 2 agents and 4 items:
  1 => {2, 3}
  2 => {1, 4}
```

## Constraints

In *constrained* fair allocation, extra criteria determine which allocations are
allowed, or *feasible*. Such constraints could, in principle, designate any
class of allocations as feasible, but the most common scenario is that the
constraints is only concerned with *bundles*, i.e., which bundles are feasible,
individually.

!!! note

    One example of a different kind of constraint is that of [Gourvès and
    Monnot](https://doi.org/10.1016/j.tcs.2018.05.018), who permit [partial
    allocations](@ref generalized), but require the set of allocated items to be
    a basis of a given matroid.

In the Allocations package, constrained instances are represented by a valuation
[`Profile`](@ref) along with a single [`Constraint`](@ref) object. (If one
wishes to combine constraints, this can be achieved using
the [`Constraints`](@ref) type.) For example, one might use a graph to indicate
items that are in conflict, and cannot be allocated to the same agent:

```jldoctest intro
julia> using Graphs

julia> V = Additive([1 1 0; 0 0 2])
Additive{Matrix{Int64}} with 2 agents and 3 items:
 1  1  0
 0  0  2

julia> G = path_graph(3)
{3, 2} undirected simple Int64 graph

julia> alloc_mnw(V, Conflicts(G)).alloc
Allocation with 2 agents and 3 items:
  1 => {2}
  2 => {1, 3}
```

In this example, agent `1` wants the first two items, while agent `2` want the
third. Allocating accordingly would provide a perfectly fair distribution, with
identical utility for both agents. However, our conflict graph—a path of length
3—tells us that no agent can have both items `1` and `2` (nor both `2` and `3`).
Our only choice, then (to avoid a Nash welfare of zero) is for the second agent
to get the first and third items.

For an overview of several types of constraints ([some of which](@ref types) are
implemented by the Allocations library), see, e.g., the survey by
[Suksompong](https://doi.org/10.1145/3505156.3505162).

### Symmetry

We call a constraint *symmetric in the agents* (or simply *symmetric*) if the
set of feasible bundles is the same for all agents. This is the most common
case, though there are exceptions. For example, [Dror et
al.](https://ojs.aaai.org/index.php/AAAI/article/view/16670) let agents have
different matroid constraints, and with so-called [*budget
constraints*](https://arxiv.org/abs/2012.03766), the agents may have different
budgets. The [`Permitted`](@ref), [`Forbidden`](@ref) and [`Required`](@ref)
constraints are also asymmetric, in general.

For many purposes, symmetry makes no difference. One important exception is
finding MMS-allocations. For symmetric constraints (and in the unconstrained
case), an agent's *maximin share* (MMS) is the greatest value she can get if she
gets the worst bundle in a feasible partition of the items, and an MMS
allocation is one where each agent gets her MMS.

If the constraint is asymmetric, however, this definition of MMS might not make
sense. The idea is that the agent partitions the items, and then gets to choose
last—i.e., gets the worst bundle, according to her. However, different
constraints apply to the different bundles, she could find a partition where she
isn't *permitted* the worst bundle, for example. There are several ways of
approaching this situation (e.g., enforcing the agent's constraint *in addition*
to the other agents' bundles, or making sure the worst bundle is feasible for
the agent in question).

Our (default) approach is to assume that the partitioning agent *doesn't know
which bundles the others would choose*. She doesn't know their constraints, and
so cannot apply them, and she doesn't know which bundle will be left for her, so
she must apply her constraint to every bundle. In other words, finding the MMS
becomes a maximin allocation problem where both the valuation function and the
constraint of the agent in question is replicated, creating $n$ "clone agents".
(For symmetric constraints, this is, of course, equivalent to the standard
approach.)

The discrepancy in constraints between finding the MMS of each agent and finding
an MMS allocation means that we can have situations where feasible allocations
exist, but the MMS of one or more agents is not defined (because no feasible MMS
partition exists). Because the duplicated constraints here represent the agent's
ability to receive any of the bundles in the partition, and not the actual
feasibility of the partition as an allocation, one strategy is to permit
leaving some items out (cf. [the next section](@ref generalized)) when finding
the MMS partitions. See [`alloc_mms`](@ref) for more on this.

## [Generalized allocations](@id generalized)

Ordinarily, an allocation is assumed to be a partition of the items being
allocated, with each item belonging to exactly one agent. `Allocation` objects,
however, are capable of representing more general forms of allocations,
where each item can belong to an arbitrary subset of the agents:

```jldoctest intro
julia> A = Allocation(1 => 1, 2 => 1, 2 => 3)
Allocation with 2 agents and 3 items, 1 unallocated:
  1 => {1}
  2 => {1, 3}
```

For [methods relying on a mixed-integer programming solver](@ref MIP-based),
such generalized allocations may be permitted by using the keyword arguments
`min_bundle`, `max_bundle`, `min_owners` and `max_owners` (as described in the
documentation of each method).

Allocating items to multiple agents is relevant, for example, in scenarios like
reviewing papers for an academic conference (as discussed by [Lesca and
Perny](https://doi.org/10.3233/978-1-60750-606-5-393)), or, perhaps, allocating
shifts to workers. Such scenarios can generally also be represented by
duplicating the items, and adding a [cardinality constraint](@ref Counts),
permitting each agent at most one copy of any item.

Allowing some items to remain unallocated (often thought of as being given to
charity) can make a big difference to the problem being solved. Some
[constraints](@ref Constraint) may make it impossible to allocate all the items.
For example, if a [conflict graph](@ref Conflicts) with a high maximum degree is
imposed, it may simply be impossible to allocate all items (cf. [Hummel and
Hetland](https://arxiv.org/abs/2104.06280)). If the feasible bundles of any
agent form an [independence
system](https://en.wikipedia.org/wiki/Independence_system) (i.e., an empty
bundle is feasible, and adding items can never repair an infeasible bundle),
permitting charity ensures that an allocation is always possible.[^2]

[^2]: The same is true of envy-based measures such as
    [EFX](https://en.wikipedia.org/wiki/Envy-free_item_allocation), which may not always
    exist even for unconstrained instances. (This is an open question at the
    time of writing.)

When allocating goods with fairness criteria that maximize efficiency in some
way (such as [leximin](@ref alloc_lmm), [MNW](@ref alloc_mnw) or [GGI](@ref
alloc_ggi)), there will be a tendency towards not needlessly leaving items to
charity. For purely envy-based criteria, however (such as [EF1](@ref alloc_ef1)
or [EFX](@ref alloc_efx)), one might in principle end up satisfying the
criterion by simply not allocating anything. In these cases, it might be better
to combine with an efficiency optimization (such as, e.g.,
[`alloc_mnw_ef1`](@ref)).

[Maximin](@ref alloc_mm) and [MMS](@ref alloc_mms) are in a special position,
here, in that they may leave items unallocated unnecessarily (i.e., not forced
by a constraint), but simply distributing these items arbitrarily will not
affect the validity of the solution. If, for some reason, charity were permitted
with any of these in an unconstrained setting, any remaining items could simply
be distributed using [`fill_random!`](@ref) or [`fill_even!`](@ref).

!!! note

    If you use [maximin](@ref alloc_mm) or [MMS](@ref alloc_mms) with a
    [constraint](@ref Constraint), distributing any leftover items is generally
    intractable, and probably best handled by the same solver that found the
    allocation to begin with. For maximin, you might consider simply switching
    to [leximin](@ref alloc_lmm) (which will also return a maximin
    allocation).[^3]

[^3]: For MMS, there is currently no ready-made solution to ensure efficiency,
    though one could, for example, perform a two-step optimization (along the
    lines of `lex_optimize!`, found in `util.jl`). First, one would use
    `alloc_mms` and `mms_alpha` to find the MMS values of all agents and the
    proportion of the MMS guarantee attainable, and then one would add the
    appropriate constraints to ensure this proportion while maximizing some
    efficiency measure, e.g., the sum of utilities.

## Solver configuration

Several allocation functions use mixed-integer linear programming via
[JuMP](https://jump.dev). Depending on the choice of MIP solver, solving even
moderately-sized instances may take a significant amount of time. Choosing a
different solver (from the default `HiGHS.Optimizer`) may speed things up
considerably. For example, with the appropriate license, one could use use
[Gurobi](https://www.gurobi.com) as follows:[^4]

[^4]: If you're a student or a researcher, Gurobi is [available for free under an academic license](https://www.gurobi.com/academia).

```julia
Allocations.conf.MIP_SOLVER = Gurobi.Optimizer
```

It is also possible to supply the `Optimizer` (or other optimizer factories,
e.g., constructed using `optimizer_with_attributes`) as the `solver` keyword
argument to the relevant allocation functions.

Normally, the MIP solvers will print out quite a lot of information about what
they're doing. If you're not interested in this output, you can generally turn
it off using some solver-specific flag, supplied to
`optimizer_with_attributes`.[^5] This is also where you'd supply other
parameters, e.g., indicating time limits, acceptable inaccuracies, etc. For
example:[^6]

[^5]: There is also the `JuMP.set_silent` function, but it requires access to the MIP model.

[^6]: See the [Gurobi](https://www.gurobi.com/documentation/10.0/refman/parameters.html) manual for explanations.

```julia
Allocations.conf.MIP_SOLVER = optimizer_with_attributes(
    Gurobi.Optimizer,
    "LogToConsole" => 0,     # No console output
    "TimeLimit" => 60,       # Finish within 60 seconds
    "MipGap" => 0.05,        # Permit 5% suboptimality
)
```

If you're unable to get rid of the output using solver parameters, a simple
solution is to just silence all output while allocating:

```julia-repl
julia> redirect_stdout(devnull) do
           alloc_mnw(V)
       end
```

If that doesn't do the trick, you could add `redirect_stderr` as well.
