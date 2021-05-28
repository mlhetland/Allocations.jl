import Base: firstindex, getindex, in, iterate, lastindex, length, show,
    summary, show_vector

using DataStructures


##############################################################################


"""
    na(X)

The number of agents represented by (e.g., valuation or allocation) `X`.
"""
function na end


"""
    agents(X)

Returns the set of agents associated with (e.g., valuation or allocation)
`X`), as an iterable of `Int`s.
"""
agents(X) = 1:na(X)


"""
    ni(X)

The number of items represented by (e.g., valuation or allocation) `X`.
"""
function ni end


"""
    items(V::Valuation)

Returns the set of items associated with (e.g., valuation or allocation) `X`,
as an iterable of `Int`s.
"""
items(X) = 1:ni(X)


# Used by Allocation and Additive
function show_agents_and_items(io::IO, X)
    n, m = na(X), ni(X)
    print(io,
        n, (n == 1 ? " agent" : " agents"),
        " and ",
        m, (m == 1 ? " item" : " items"))
end


## Allocations ###############################################################


"""
    struct Allocation <: Any

A mapping `A` from agents `i` to their assigned bundles `bundle(A, i)`. Agents
and items are represented as `Int`s, and bundles as `Set`s of `Int`s. The
`Allocation` also maintains an inverse mapping, from items `g` to their set of
owners, `owners(A, g)`. To keep these in sync, the bundles should be modified
using `give!` and `deny!`, rather than altering the bundle sets directly.
"""
struct Allocation
    # (Might make more sense to use `bundles` here, but to be consistent, we'd
    # then also need to rename `owners`...)
    bundle::Vector{Set{Int}}    # bundle[i]: Agent i's bundle
    owners::Vector{Set{Int}}    # owners[g]: Item g's owners
end


"""
    Allocation(n::Int, m::Int)

Construct an empty allocation with `n` agents and `m` items.
"""
Allocation(n::Int, m::Int) =
    Allocation([Set{Int}() for i = 1:n], [Set{Int}() for i = 1:m])


bundle(A) = A.bundle
owners(A) = A.owners


na(A::Allocation) = length(bundle(A))
ni(A::Allocation) = length(owners(A))


## Allocation printing


# Rather than reinventing printing code, we'll let built-in functions for
# AbstractDict and AbstractVector much of the heavy lifting, using a couple of
# one-off wrappers (AllocShowWrap and BundleShowWrap).


function show_bundle(io::IO, S)
    seq = sort(collect(S))
    # Supplying :typeinfo to "pretend" we've already printed the type info
    ctx = IOContext(io, :typeinfo => typeof(seq))
    show_vector(ctx, seq, "{", "}")
end


# Single-use -- assumes alloc will not be modified after construction, because
# this won't be reflected in pairs. This wrap is used to leverage the
# MIME"text/plain" printing of an AbstractDict; it is not used with single-line
# printing of Allocations (which is implemented more directly).
struct AllocShowWrap <: AbstractDict{Int,Set{Int}}
    alloc::Allocation
    pairs
end
AllocShowWrap(A) =
    AllocShowWrap(A, enumerate(BundleShowWrap.(A.bundle)))


summary(io::IO, w::AllocShowWrap) = summary(io, w.alloc)
show(io::IO, w::AllocShowWrap) = show(io, w.alloc)
length(w::AllocShowWrap) = na(w.alloc)
iterate(w::AllocShowWrap, args...) = iterate(w.pairs, args...)


struct BundleShowWrap
    bundle
end


show(io::IO, w::BundleShowWrap) = show_bundle(io, w.bundle)


function summary(io::IO, A::Allocation)
    print(io, "Allocation with ")
    show_agents_and_items(io, A)
    u = length([g for g in items(A) if !owned(A, g)])
    if u != 0
        print(io, ", ", u, " unallocated")
    end
end


function show(io::IO, A::Allocation)
    seq = BundleShowWrap.(A.bundle)
    # Supplying :typeinfo to "pretend" we've already printed the type info
    ctx = IOContext(io, :typeinfo => typeof(seq))
    show_vector(ctx, seq)
end


show(io::IO, m::MIME"text/plain", A::Allocation) = show(io, m, AllocShowWrap(A))


## Allocation accessors and manipulators


"""
    bundle(A, i)

The set of items allocated to agent `i` in the allocation `A`. The returned
`Set` should be treated as read-only.
"""
bundle(A, i) = bundle(A)[i]


"""
    owners(A, g)

The set of agents to which item `g` has been allocated in the allocation `A`.
The returned `Set` should be treated as read-only.
"""
owners(A, g) = owners(A)[g]


"""
    owner(A, g)

The agent to which item `g` has been allocated in the allocation `A`. Will
produce an error if `g` has been allocated to more than one agent.
"""
owner(A, g) = only(owners(A, g))


"""
    owned(A, g)

Whether or not the item `g` is owned by any agent in the allocation `A`.
"""
owned(A, g) = !isempty(owners(A, g))


"""
    give!(A, i, g::Int)

Give agent `i` the object `g` in the `Allocation` `A`.
"""
function give!(A, i, g::Int)
    push!(bundle(A, i), g)
    push!(owners(A, g), i)
    return A
end


"""
    give!(A, i, B)

Give agent `i` the bundle `B` in the `Allocation` `A`.
"""
function give!(A, i, B)
    union!(bundle(A, i), B)
    for g in B
        push!(owners(A, g), i)
    end
    return A
end


"""
    deny!(A, i, g)

Deny agent `i` the object `g`, which it has previously been given, in the
allocation `A`.
"""
function deny!(A, i, g)
    delete!(bundle(A, i), g)
    delete!(owners(A, g), i)
    return A
end


"""
    fill_even!(A)

Fill out the allocation by distributing the unallocated items evenly, by
repeatedly giving the next unallocated item to the agent with the fewest items
(ties broken arbitrarily).
"""
function fill_even!(A)
    n, m = na(A), ni(A)
    pq = PriorityQueue{Int, Int}()
    for i = 1:n
        enqueue!(pq, i, length(bundle(A, i)))
    end
    for g = 1:m
        owned(A, g) && continue
        i, u = peek(pq)
        give!(A, i, g)
        pq[i] = u + 1
    end
    return A
end


## Valuations ################################################################


"""
    abstract struct Valuation <: Any

An abstract type representing a valuation oracle. Which functions are used to
query it depends on the kind of valuation functions it represents. Additive
valuations act on individual objects, and simply sum those values over a
bundle, but oracles with quite different kinds of queries are possible for
valuations with other properties (see, e.g., [Fair Allocation of Indivisible
Goods: Improvements and
Generalizations](https://dl.acm.org/doi/10.1145/3219166.3219238) by Ghodsi et
al., 2018).
"""
abstract type Valuation end


"""
     Allocation(V::Valuation)

Construct an empty allocation with a number of agents and items equal to that of
the instance `V`.
"""
Allocation(V::Valuation) = Allocation(na(V), ni(V))


"""
    value(V::Valuation, i, S)
    value(V::Valuation, i, g::Int)

The value agent `i` places on bundle `S`, according to the oracle `V`. The
second form is a shortcut for `value(V, i, [g])`, which will generally be more
efficient. Note that the value of `S` may *not* in general be the sum of the
values of its items; that property is unique to `Additive` valuations.
"""
function value end


"""
    value(V::Valuation, i, A::Allocation)

The value agent `i` receives in allocation `A`, under the valuation `V`.
"""
value(V::Valuation, i, A::Allocation) = bundle_value(V, i, A)

# Use for disambiguation definitions for subtypes:
bundle_value(V, i, A) = value(V, i, bundle(A, i))


"""
    value_1(V::Valuation, i, S)

The value agent `i` places on bundle `S`, *up to one item*, that is, the
smallest value `i` can place on bundle `S` after removing (at most) one item,
according to the oracle `V`.
"""
function value_1 end


"""
    value_x(V::Valuation, i, S)

The value agent `i` places on bundle `S`, *up to any item*, that is, the
largest value `i` can place on bundle `S` after removing one item (or no
items, if the bundle is empty), according to the oracle `V`.
"""
function value_x end


"""
    isintegral(V::Valuation)

Test whether every value provided by `V` is an integer.
"""
function isintegral end


"""
    isnonnegative(V::Valuation)

Test whether every value provided by `V` is nonnegative.
"""
function isnonnegative end


"""
    matrix(V::Valuation)

Return a matrix `X` where `X[i, g]` is `value(V, i, g)`. May not be very useful
in general (especially if calculating single-item values isn't efficient to
begin with), but if such a matrix is available as part of the valuation
implementation (as with `Additive`), it may be returned directly.
"""
matrix(V::Valuation) = [value(V, i, g) for i in agents(V), g in items(V)]


"""
    struct Additive{T <: AbstractMatrix} <: Valuation

An additive valuation oracle, representing how each agent values all possible
bundles. Because of additivity, this is easily "lifted" from the values of
individual items, by addition, with an empty bundle given a value of zero. By
default, the valuation is constructed from a real matrix `X`, supplied to the
default constructor, where `X[i, g]` is agent `i`'s value for item `g`.
"""
struct Additive{T <: AbstractMatrix{<:Real}} <: Valuation
    values::T
end


"""
    Additive(n, m)

Create an additive valuation for `n` agents and `m` items where all values are
set to zero.
"""
Additive(n, m) = Additive(zeros(n, m))


"""
    Valuation(X::Matrix)

Alias for `Additive(X)`.
"""
Valuation(X::Matrix) = Additive(X)


na(V::Additive) = size(V.values, 1)
ni(V::Additive) = size(V.values, 2)


function summary(io::IO, A::Additive{T}) where {T}
    print(io, "Additive{$T}")
end


function show(io::IO, ::MIME"text/plain", A::Additive)
    summary(io, A)
    print(io, " with ")
    show_agents_and_items(io, A)
    println(":")
    Base.print_matrix(io, A.values)
end


isintegral(V::Additive) = all(isinteger.(V.values))


isnonnegative(V::Additive) = all(V.values .>= zero(eltype(V.values)))


"""
    value(V::Additive, i, g::Int)

The value of item `g`, according to agent `i`.
"""
value(V::Additive, i, g::Int) = V.values[i, g]


# Disambiguation:
value(V::Additive, i, A::Allocation) = bundle_value(V, i, A)


# The bundle value is "lifted" from item values by addition.
# TODO When Julia 1.6 comes out, can use init keyword argument
value(V::Additive, i, S) =
    isempty(S) ? zero(eltype(V.values)) : sum(value(V, i, g) for g in S)


# For the additive case, we can just subtract the highest item value.
# TODO When Julia 1.6 comes out, can use init keyword argument
value_1(V::Additive, i, S) =
    value(V, i, S) -
    (isempty(S) ? zero(eltype(V.values)) : maximum(value(V, i, g) for g in S))


# For the additive case, we can just subtract the lowest item value.
# TODO When Julia 1.6 comes out, can use init keyword argument
value_x(V::Additive, i, S) =
    value(V, i, S) -
    (isempty(S) ? zero(eltype(V.values)) : minimum(value(V, i, g) for g in S))


"""
    value!(V::Additive, i, g::Int, v)

Set the value of item `g`, according to agent `i`, to `v`.
"""
value!(V::Additive, i, g, v) = V.values[i, g] = v


"""
    matrix(V::Additive)

Return the underlying valuation matrix of `V`.
"""
matrix(V::Additive) = V.values


"""
    normalize(V)

Scale the values of `V` such that ``v_i(M) = n`` for all agents ``i``.
"""
normalize(V::Additive) =
    Additive(V.values .* [na(V) / value(V, i, items(V)) for i in agents(V)])


## Constraints ###############################################################


"""
    abstract type Constraint <: Any

Abstract supertype of various kinds of constraints. An allocation problem is
assumed to consist of a `Valuation` object and at most one `Constraint`
object, embodying any and all constraints placed on feasible solutions.
"""
abstract type Constraint end


"""
    mutable struct Category

One of the categories in a `Counts` constraint, from which each agent can hold
at most a given number of items. The category supports iteration (over its
members), and the threashold is available through the `threshold` accessor.
"""
mutable struct Category
    members::Set{Int}
    threshold::Int
end


iterate(c::Category, args...) = iterate(c.members, args...)
length(c::Category) = length(c.members)


"""
    threshold(c::Category)

The maximum number of items any agent can receive from the given category, as
part of a `Counts` constraint.
"""
threshold(c::Category) = c.threshold


"""
    mutable struct OrderedCategory

Used in place of `Category` when handling an ordered instance. The instance is
assumed to be such that items in the range `index:index + n_items - 1` belong to
the given category, i.e., the items of a category occupy a contiguous range of
integers.
"""
mutable struct OrderedCategory
    index::Int
    n_items::Int
    threshold::Int
end


in(g, c::OrderedCategory) = c.index <= g < c.index + c.n_items
lastindex(c::OrderedCategory) = length(c)
length(c::OrderedCategory) = c.n_items

getindex(c::OrderedCategory, i::Int) =
    getindex(c.index:c.index + c.n_items - 1, i)

getindex(c::OrderedCategory, v::UnitRange{Int64}) =
    getindex(c.index:c.index + c.n_items - 1, v)

iterate(c::OrderedCategory, args...) =
    iterate(c.index:c.index + c.n_items - 1, args...)


"""
    floor_n(c::OrderedCategory, n)

One `n`th of the number of items in the category, rounded down.
"""
floor_n(c::OrderedCategory, n) = floor(Int, c.n_items/n)


"""
    ceil_n(c::OrderedCategory, n)

One `n`th of the number of items in the category, rounded up.
"""
ceil_n(c::OrderedCategory, n) = ceil(Int, c.n_items/n)


"""
    required(c::OrderedCategory, n)

The number of items the next agent must take in order to keep the instance
valid, i.e., for there to be a maximum of `(n - 1) * threshold` remaining
items.
"""
required(c::OrderedCategory, n) = max(c.n_items - (n - 1) * c.threshold, 0)


"""
    struct Counts{T} <: Constraint

The *cardinality constraints* introduced by Biswas and Barman in their 2018
paper [Fair Division Under Cardinality
Constraints](https://www.ijcai.org/proceedings/2018/13). This is a form of
constraint consisting of several `Category` objects, available through
indexing or iteration. Any agent may hold at most a given number of items from
any given category.
"""
struct Counts{T} <: Constraint
    categories::Vector{T}
end


"""
    Counts(args::Pair...)

Create a `Counts` object where each pair `x => k` becomes a category with
members `Set(x)` and threshold `k`.
"""
Counts(args::Pair...) = Counts([Category(Set(p[1]), p[2]) for p in args])


getindex(C::Counts, i) = C.categories[i]
length(C::Counts) = length(C.categories)
iterate(C::Counts, args...) = iterate(C.categories, args...)


"""
    struct Conflicts{T <: AbstractGraph} <: Constraint

A kind of constraint -- or set of constraints -- that indicates that certain
items conflict, and thus cannot be allocated to the same agent. The
constraints are represented as a *conflict graph*
(`LightGraphs.AbstractGraph`), with items as nodes, and edges representing
conflicts. The `Conflicts` type is just a wrapper for dispatch purposes, with
the underlying graph available through the `graph` accessor.
"""
struct Conflicts{T <: AbstractGraph} <: Constraint
    graph::T
end


"""
    graph(C::Conflicts)

Return the conflict graph wrapped by a `Conflicts` object.
"""
graph(C::Conflicts) = C.graph


## Reductions ################################################################


"""
    mutable struct Reduction{S, T, U, V}

A reduction from one instance of a fair allocation problem to another. Contains
information about the valuations in the reduced instance, through an object of
type `S`. There must exist functions `agents(s::S)` and `items(s::S)` that
return an iterator of, respectively, the agents and items in the reduced
instance. The reduction can also contain information about the constraints in
the reduced instance, through an object of type `T`.

In addition, the reduction contains two mappings, `λi` (of type `U`) and `λg`
(of type `V`). Both types should be indexable (for `i ∈ agents(s)` and `g ∈
items(s)`, respectively). `λi[i]` and `λg[g]` should return the agent and item
identifier in the original instance of, respectively, agent `i` and item `g` in
the reduced instance.

The reduction also contains a function that can convert an allocation in the
reduced instance to one in the original instance.

The default constructor is `Reduction(V, C, λi, λg, transform::Function)`.
"""
mutable struct Reduction{S, T, U, V}
    V::S
    C::T
    λi::U
    λg::V
    transform::Function
end


"""
    Reduction(V, λi, λg, transform)

A simplified constructor for when there are no constraints.
"""
Reduction(V, λi, λg, transform) = Reduction(V, nothing, λi, λg, transform)


"""
    Reduction(V, C)

A simplified constructor for when either no changes have been performed or
changes only concern the valuations and/or constraints.
"""
Reduction(V, C) = Reduction(V, C, agents(V), items(V), identity)


"""
    Reduction(V)

A simplified constructor for when either no changes have been performed or
changes only concern the valuations.
"""
Reduction(V) = Reduction(V, nothing)


"""
    valuation(R::Reduction)

Returns the valuations in the reduced instance.
"""
valuations(R::Reduction) = R.V


"""
    constraints(R::Reduction)

Returns the constraints in the reduced instance
"""
constraints(R::Reduction) = R.C


"""
    transform(R::Reduction, A::Allocation)

Converts the given allocation for the reduced instance to one for original
instance. The way the convertion occurs depends on the given reduction.
"""
transform(R::Reduction, A::Allocation) = R.transform(A)


"""
    agent(R::Reduction, i)

Converts the agent identifier `i` from the reduced instance to the agent
identifier of the same agent in the original instance.
"""
agent(R::Reduction, i) = R.λi[i]


"""
    item(R::Reduction, g)

Converts the item identifier `g` from the reduced instance to the item
identifier of the same item in the original instance.
"""
item(R::Reduction, g) = R.λg[g]


"""
    chain(R₁::Reduction, R₂::Reduction)

Assumes that R₂ is a reduction of the reduced instance of R₁. Combines the two
reductions, so that the original instance is the original instance of R₁ and the
reduced instance is the reduced instance of R₂ (essentially diagram-order
composition of the reductions).
"""
function chain(R₁::Reduction, R₂::Reduction)
    return Reduction(
            R₂.V, R₂.C,
            [agent(R₁, agent(R₂, i)) for i in agents(R₂.V)],
            [item(R₁, item(R₂, g)) for g in items(R₂.V)],
            (A) -> transform(R₁, transform(R₂, A))
        )
end
