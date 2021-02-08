import Base: getindex, iterate, length, show, summary


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


function show_agents_and_items(io::IO, X)
    n, m = na(X), ni(X)
    print(io,
        n, (n == 1 ? " agent" : " agents"),
        " and ",
        m, (m == 1 ? " item" : " items"))
end


show_bundle(io::IO, S) =
    print(io, "{", join(sort(collect(S)), ", "), "}")


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


function summary(io::IO, A::Allocation)
    print(io, "Allocation with ")
    show_agents_and_items(io, A)
    u = length([g for g in items(A) if isempty(owners(A, g))])
    if u != 0
        print(io, ", ", u, " unallocated")
    end
end


function show(io::IO, A::Allocation)
    print(io, "[",
        join([sprint(show_bundle, bundle(A, i)) for i in agents(A)], ", "),
    "]")
end


function show(io::IO, ::MIME"text/plain", A::Allocation)
    summary(io, A)
    print(io, ":")
    for i in agents(A)
        print(io, "\n  ", i, " => ")
        show_bundle(io, bundle(A, i))
    end
end


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
    give!(A, i, g)

Give agent `i` the object `g` in the `Allocation` `A`.
"""
function give!(A, i, g)
    push!(bundle(A, i), g)
    push!(owners(A, g), i)
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


## Valuations ################################################################


"""
    struct Valuation <: Any

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
    value(V::Valuation, i, S)

The value agent `i` places on bundle `S`, according to the oracle `V`.
"""
function value end


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
    value!(V::Additive, i, j::Int, v)

Set the value of item `j`, according to agent `i`, to `v`.
"""
value!(V::Additive, i, j, v) = V.values[i, j] = v


## Constraints ###############################################################


"""
    abstract type Constraint <: Any

Abstract supertype of various kinds of constraints. An allocation problem is
assumed to consist of a `Valuation` object and at most one `Constraint`
object, embodying any and all constraints placed on feasible solutions.
"""
abstract type Constraint end


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
