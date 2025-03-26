abstract type Matroid end


Base.show(io::IO, M::Matroid) =
    print(io, "$(nameof(typeof(M)))($(M.n), $(rank(M)))")

Base.show(io::IO, ::MIME"text/plain", M::Matroid) =
    print(io, "$(nameof(typeof(M))) with size $(M.n) and rank $(rank(M))")


struct ClosedSetsMatroid <: Matroid
    n::Integer # Size of universe
    r::Integer # Final rank (r == length(F)).
    F::Vector{Set{SmallBitSet}} # Closed sets by rank
end

Base.copy(M::ClosedSetsMatroid) = ClosedSetsMatroid(M.n, M.r, M.F)


struct FullMatroid <: Matroid
    n::Integer
    r::Integer
    F::Vector{Set{SmallBitSet}} # Closed sets by rank
    I::Vector{Set{SmallBitSet}} # Independent sets by rank
    C::Union{Nothing,Set{SmallBitSet}} # Circuits
end

Base.copy(M::FullMatroid) = FullMatroid(M.n, M.r, M.F, M.I, M.C)


struct GraphicMatroid <: Matroid
    g::Graph
    n::Integer
    r::Integer
end

GraphicMatroid(g::Graph) = GraphicMatroid(g, ne(g), length(kruskal_mst(g)))
Base.copy(M::GraphicMatroid) = GraphicMatroid(M.g, M.n, M.r)


struct UniformMatroid <: Matroid
    n::Integer
    r::Integer
end

Base.copy(M::UniformMatroid) = UniformMatroid(M.n, M.r)


struct PartitionMatroid <: Matroid
    n::Integer
    categories::Vector{UnitRange{Integer}}
end

Base.copy(M::PartitionMatroid) = PartitionMatroid(M.n, M.categories)


FreeMatroid(n) = UniformMatroid(n, n)
ZeroMatroid(n) = UniformMatroid(n, 0)


ground_set(M::ClosedSetsMatroid) = SmallBitSet(1:M.n)
ground_set(M::FullMatroid) = SmallBitSet(1:M.n)
ground_set(M::UniformMatroid) = BitSet(1:M.n)
ground_set(M::GraphicMatroid) = BitSet(1:M.n)


"""
    is_indep(M::Matroid, S)

Determine whether the set `S` is independent in the matroid `M`.
"""
function is_indep end


# Determines whether a given set S is independent in the matroid M, given by the
# closed sets of M grouped by rank. Uses (I.1) in Greene (1991).
function is_indep(M::ClosedSetsMatroid, S)
    t = length(S)

    if t == 0
        return true
    elseif t > M.r
        return false
    end

    for F in M.F[t]
        if issubset(S, F)
            return false
        end
    end

    return true
end


function is_indep(M::FullMatroid, S)
    card = length(S)
    if card == 0
        return true
    elseif card + 1 > length(M.I)
        return false
    end
    return S in M.I[card+1]
end


is_indep(M::UniformMatroid, S) = length(S) <= M.r


function is_indep(M::GraphicMatroid, S)
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    subgraph, _vmap = induced_subgraph(M.g, edgelist)
    return !is_cyclic(subgraph)
end


"""
    rank(M::Matroid, S)

Returns the rank of the set `S` in `M`, i.e., the size of the largest
independent subset of `S`.
"""
function rank end


rank(M::Matroid) = M.r


function rank(M::ClosedSetsMatroid, S)
    for (r, Fr) in enumerate(M.F), B ∈ Fr
        if issubset(S, B)
            return r - 1
        end
    end
end


function rank(M::FullMatroid, S)
    for (r, Fr) in enumerate(M.F), F ∈ Fr
        if issubset(S, F)
            return r - 1
        end
    end
end


rank(M::UniformMatroid, S) = min(length(S), M.r)


rank(M::GraphicMatroid, S) = length(S) > 0 ? length(minimal_spanning_subset(M, S)) : 0


"""
    is_circuit(M::Matroid, S)

Determines whether `S` is a circuit in `M` (i.e., `rank(M, S) == length(S) - 1`).
"""
function is_circuit end


# Determines whether `S` is a circuit in the matroid, using (C.1) and (C.2) in
# Greene (1991).
function is_circuit(M::ClosedSetsMatroid, S)
    t = length(S)

    # (C.1) S ⊆ F for some F ∈ F_{t-1}.
    # t equals t-1 due to 1-indexing.
    any(F -> issubset(S, F), M.F[t]) || return false

    # (C.2) |S ∩ F| ≤ r(F) for all F ∈ F_{t-2}.
    # t-1 equals t-2 due to 1-indexing.
    return all(F -> length(intersect(S, F)) <= t - 2, M.F[t-1])
end


is_circuit(M::FullMatroid, S) = S in M.C


is_circuit(M::UniformMatroid, S) = M.r == length(S) - 1


function is_circuit(M::GraphicMatroid, S)
    return rank(M, S) == length(S) - 1
end


"""
    minimal_spanning_subset(M::Matroid, S)

Finds a minimal spanning subset of `S` in `M`. If `S` is the ground set of `M`,
this produces a basis of `M`.
"""
function minimal_spanning_subset end


minimal_spanning_subset(M::ClosedSetsMatroid, S) = _mss(M, 0, S)
minimal_spanning_subset(M::FullMatroid, S) = _mss(M, 0, S)


# Algorithm 3.1 from Greene (1991)
function _mss(M::Union{ClosedSetsMatroid,FullMatroid}, j::Integer, A)
    B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]

    while length(B) == 0
        if j >= length(A) - 1
            return A
        end

        j += 1
        B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]
    end

    _mss(M, j, setdiff(A, rand(reduce(union, B))))
end


minimal_spanning_subset(M::UniformMatroid, S) = Set(collect(S)[1:min(end, M.r)])


# Uses Kruskal's algorithm to find a minimal spanning tree over M.G.
function minimal_spanning_subset(M::GraphicMatroid, S)
    # Maintain a map for converting edges to the corresponding items, while
    # ignoring the order of the vertices in the edge
    edgemap = Dict(Set([e.src, e.dst]) => i
                   for (i, e) in enumerate(edges(M.g)))
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    subgraph, vmap = induced_subgraph(M.g, edgelist)
    mst = kruskal_mst(subgraph)
    # Convert edges back into items
    return Set([edgemap[Set([vmap[e.src], vmap[e.dst]])] for e in mst])
end


"""
    minimal_spanning_subsets(M::Matroid, S)

Finds all minimal spanning subsets of `S` in `M`. If `S` is the ground set of
`M`, this finds all the bases of `M`.

See also [`bases`](@ref).
"""
function minimal_spanning_subsets end


minimal_spanning_subsets(M::ClosedSetsMatroid, S) = _mss_all(M, 0, S)
minimal_spanning_subsets(M::FullMatroid, S) = _mss_all(M, 0, S)


# A modification of Algorithm 3.1 from Greene (1991) that finds all minimal
# spanning subsets of A ⊆ E.
function _mss_all(M, j::Integer, A)
    B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]

    while length(B) == 0
        if j >= length(A) - 1
            # Make sure to add A to a new set, without flattening (since A is
            # already a set)
            return Set([A])
        end

        j += 1
        B = [intersect(F, A) for F in M.F[j+1] if length(intersect(F, A)) > j]
    end

    bases = Set{AbstractSet{Int}}()
    T = reduce(union, B)
    while !isempty(T)
        x = minimum(T)
        union!(bases, _mss_all(M, j, setdiff(A, x)))
        setdiff!(T, x)
    end

    return bases
end


function minimal_spanning_subsets(M::UniformMatroid, S)
    M.n < 64 || throw(DomainError(M.n, "matroid size >= 64 is not supported"))
    len = min(length(S), M.r)
    mask = set_to_bits(S)
    return Set([bits_to_set(i) for i in (2^len - 1):(2^M.n - 1)
                if count_ones(i) == len && i & mask == i])
end


minimal_spanning_subsets(M::GraphicMatroid, S) = throw("unimplemented")


"""
    bases(M::Matroid)

Finds the set of bases of the matroid `M`.
"""
function bases end


bases(M::ClosedSetsMatroid) = _mss_all(M, 0, ground_set(M))
bases(M::FullMatroid) = _mss_all(M, 0, ground_set(M))


function bases(M::UniformMatroid)
    M.n < 64 || throw(DomainError(M.n, "matroid size >= 64 is not supported"))
    return Set([bits_to_set(i) for i in 1:(2^M.n - 1) if count_ones(i) == M.r])
end


bases(M::GraphicMatroid) = throw("unimplemented")


"""
    closure(M::Matroid, S)

Returns the closure of a set `S` in the matroid `M`.

Given a set `S ⊆ E`, add all elements `x ∈ E` such that `S` has no increase in
rank.
"""
function closure end


function closure(M::ClosedSetsMatroid, S)
    for Fr in M.F, B in Fr
        if issubset(S, B)
            return B
        end
    end
end


function closure(M::FullMatroid, S)
    for Fr in M.F, B in Fr
        if issubset(S, B)
            return B
        end
    end
end


closure(M::UniformMatroid, S) = length(S) < M.r ? S : ground_set(M)


function closure(M::GraphicMatroid, S)
    edgelist = [e for (i, e) in enumerate(edges(M.g)) if i in S]
    _, vmap = induced_subgraph(M.g, edgelist)
    return Set([i for (i, e) in enumerate(edges(M.g))
                if [e.src, e.dst] ⊆ vmap || e.src == e.dst])
end


"""
    is_closed(M::Matroid, S)

Returns `true` if the given set `S` is closed in the matroid `M`, i.e., if the
closure of `S` is equal to itself.
"""
is_closed(M::Matroid, S) = issetequal(closure(M, S), S)


## Tests for the axioms for the closed sets of a matroid, as given by Knuth
## (1975).


"""
    matroid_c1(M::Union{ClosedSetsMatroid, FullMatroid})

The first axiom for a matroid defined by closed sets, as described in [Knuth's
1975 paper](https://doi.org/10.1016/0012-365X(75)90075-8):

> The ground set is closed. \$E ∈ F\$.
"""
function matroid_c1(M::Union{ClosedSetsMatroid,FullMatroid})
    last(M.F) == Set((ground_set(M),))
end


"""
    matroid_c2(M::Union{ClosedSetsMatroid, FullMatroid})

The second axiom for a matroid defined by closed sets, as described in [Knuth's
1975 paper](https://doi.org/10.1016/0012-365X(75)90075-8):

> The intersection of two closed sets is a closed set. If \$A, B ∈ F\$, then
> \$A ∩ B ∈ F\$.
"""
function matroid_c2(M::Union{ClosedSetsMatroid,FullMatroid})
    F = reduce(∪, M.F)
    for A ∈ F
        for B ∈ F
            if !(intersect(A, B) in F)
                return false
            end
        end
    end
    return true
end


"""
    matroid_c3(M::Union{ClosedSetsMatroid, FullMatroid})

The third axiom for a matroid defined by closed sets, as described in [Knuth's
1975 paper](https://doi.org/10.1016/0012-365X(75)90075-8):

> If \$A ∈ F\$ and \$a, b ∈ E - A\$, then \$b\$ is a member of all sets
> containing \$A ∪ {a}\$ if and only if \$a\$ is a member of all sets containing
> \$A ∪ {b}\$.
"""
function matroid_c3(M::Union{ClosedSetsMatroid,FullMatroid})
    E = ground_set(M)
    F = reduce(∪, M.F)
    delete!(F, SmallBitSet())
    for A ∈ F
        t1 = setdiff(E, A)
        while !isempty(t1)
            a = minimum(t1)
            t2 = setdiff(t1, a)
            while !isempty(t2)
                b = minimum(t2)
                ā = reduce(intersect, [B for B in F if issubset(union(A, a), B)])
                b̄ = reduce(intersect, [B for B in F if issubset(union(A, b), B)])

                if issubset(b, ā) != issubset(a, b̄)
                    return false
                end

                delete!(t2, b)
            end
            delete!(t1, a)
        end
    end
    return true
end
