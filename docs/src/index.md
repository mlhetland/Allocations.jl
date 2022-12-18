# Allocations.jl

The Allocations package deals with the fair allocation of indivisible items
to a set of agents. An instance of this problem consists of:

- A set $N$ of $n$ *agents* and a set $M$ of $m$ *items*;
- One *valuation function* $v_i:2^M\to\R$ for each agent $i\in N$;

For simplicity, we may simply let $N=\{1,\dots,n\}$ and $M=\{1,\dots,m\}$.[^1]

The goal is to find an *allocation* $A=(A_1,\dots,A_n)$, which gives a *bundle*
$A_i\subseteq M$ to each agent, and which satisfies certain fairness criteria.
The allocation is usually required to form a partition of $M$, though in some
scenarios, one may deviate from this.

In `Allocations.jl`, an instance is represented by the *valuation profile*
$V=\{v_i:i\in N\}$, in the form of a `Profile` object. Given a profile `V`, the
various parts of the instance may be accessed as follows:

- `items(V)`: The item set $N$, as an iterable;
- `agents(V)`: The agent set $M$, as an iterable;
- `value(V, i, S)`: The value agent i assigns to the set $S\subseteq
  M$, i.e., $v_i(S)$.

`S` may be any iterable. If it a single item (i.e., `Int`) `g`, it is
interpreted as a singleton, `[g]` (though usually handled more efficiently).
The number of agents and items, respectively, is found using `na(V)` and `ni(V)`

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

The bundles of each agent is available through the `bundle` function:

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

## Solver configuration

Several allocation functions use mixed-integer linear programming via
[JuMP](https://jump.dev). Depending on the choice of MIP solver, solving even
moderately-sized instances may take a significant amount of time. Choosing a
different solver (from the default `HiGHS.Optimizer`) may speed things up
considerably. For example, with the appropriate license, one could use use
[Gurobi](https://www.gurobi.com) as follows:[^2]

```julia
Allocations.conf.MIP_SOLVER = Gurobi.Optimizer
```

It is also possible to supply the `Optimizer` (or other optimizer factories,
e.g., constructed using `optimizer_with_attributes`) as the `solver` keyword
argument to the relevant allocation functions.

Normally, the MIP solvers will print out quite a lot of information about what
they're doing. If you're not interested in this output, you can generally turn
it off using some solver-specific flag, supplied to
`optimizer_with_attributes`.[^3] This is also where you'd supply other
parameters, e.g., indicating time limits, acceptable inaccuracies, etc. For
example:[^4]

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

[^1]: The latter is less common, presumably because the two sets then intersect. An alternative is to let $M$ be a set of opaque objects $g_j$, for $j=1\dots m$.
[^2]: If you're a student or a researcher, Gurobi is [available for free under an academic license](https://www.gurobi.com/academia).
[^3]: There is also the `JuMP.set_silent` function, but it requires access to the MIP model.
[^4]: See the [Gurobi](https://www.gurobi.com/documentation/10.0/refman/parameters.html) manual for explanations.
