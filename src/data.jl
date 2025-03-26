using Distributions
using Graphs

# Note: Although `Uniform` represents a uniform distribution on [0,1], the
# implementation (https://github.com/JuliaStats/Distributions.jl/blob/a9b0e3c9)
# actually samples from [0,1), by using rand().


"""
    rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    rand_profile(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())

Generate a random additive valuation profile (an `Additive` object) with the
number of agents and items, agents and values chosen using `rand` with the
given values as the first argument. Here `n` is used directly, while `m` should
be a *function of the number of agents*, which returns an argument for `rand`.
Similarly, `v` should be a function of both the number of agents and items.

The defaults for `n` and `m` are taken from [Hummel and
Hetland](https://arxiv.org/abs/2104.06280), who based them on based on
real-world data from [Caragiannis et al.](https://doi.org/10.1145/3355902), with
`n` at most 10, and the average `m/n` approximately 3.

The distribution for the values can be univariate, in which case it is used
independently for each value. For example, to generate Gaussian valuations, with
the [`Distributions`](https://github.com/JuliaStats/Distributions.jl) package,
use `v=(n, m)->Normal()`.

However, it can also be multivariate, in which case each sample should be a
vector of the same length as the number of items, representing the valuation
function of a single agent. For example, to get a `Dirichlet` distribution,
where an agent's values sum to `1`, you could use `v=(n, m)->Dirichlet(m, 2)`.

It is also possible to use a matrix-variate distribution, such as a matrix
normal distribution, where each sample should then be an `m` by `n` matrix, with
each column representing the valuation function of a single agent.

!!! warning

    Note that the matrix samples should be the *transpose* of the ones used in
    the resulting profile. This is to maintain consistency with the multivariate
    distributions, which produce column vectors.

`rand_profile` is an alias for `rand_additive`.
"""
function rand_additive(; n=2:10, m=n->2n:4n, v=(n, m)->Uniform(), rng=default_rng())
    nn = rand(rng, n)
    mm = rand(rng, m(nn))
    vv = v(nn, mm)
    X = rand(rng, nn, mm)
    rand!(rng, vv, X')
    return Additive(X)
end

const rand_profile = rand_additive


"""
    rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    rand_conflicts_ws98(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Watts–Strogatz model](https://doi.org/10.1038/30918). The
keyword arguments `k` and `β` specify the possible values for the corresponding
parameters \$k\$ and \$\\beta\$, which are generated using `rand`. The defaults
are taken from [Hummel and Hetland](https://arxiv.org/abs/2104.06280). Note that
the parameter \$k\$ should be an even number, which Watts and Strogatz assume to
be much smaller than \$m\$.
"""
function rand_conflicts_ws98(m; k=2:2:div(m, 2), β=Uniform(), rng=default_rng())
    G = watts_strogatz(m, rand(rng, k), rand(rng, β), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ws98(V::Profile; kwds...) = rand_conflicts_ws98(ni(V); kwds...)


"""
    rand_conflicts_er59(m, p=Uniform(), rng=default_rng())
    rand_conflicts_er59(V::Profile; ...)
    rand_conflicts(m; ...)
    rand_conflicts(m::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Erdős–Rényi
model](https://doi.org/10.5486%2FPMD.1959.6.3-4.12). The keyword argument `p`
specifies the possible values for the corresponding parameter \$p\$, which is
generated using `rand`.

`rand_conflicts` is an alias for `rand_conflicts_er59`.
"""
function rand_conflicts_er59(m; p=Uniform(), rng=default_rng())
    G = erdos_renyi(m, rand(rng, p), rng=rng)
    return Conflicts(G)
end

rand_conflicts_er59(V::Profile; kwds...) = rand_conflicts_er59(ni(V); kwds...)
const rand_conflicts = rand_conflicts_er59


"""
    rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    rand_conflicts_ba02(V::Profile; ...)

Generate a random `Conflicts` contraint, whose underlying graph is constructed
according to the [Barabási–Albert
model](https://arxiv.org/abs/cond-mat/0106096). The keyword argument `k`
specifies the possible values for the corresponding parameter \$k\$, which is
generated using `rand`.
"""
function rand_conflicts_ba02(m; k=1:m, rng=default_rng())
    G = barabasi_albert(m, rand(rng, k), rng=rng)
    return Conflicts(G)
end

rand_conflicts_ba02(V::Profile; kwds...) = rand_conflicts_ba02(ni(V); kwds...)


"""
    rand_matroid_er59(m; n=nothing, rng=default_rng())
    rand_matroid_er59(V::Profile; ...)

Generate a random `GraphicMatroid`, whose underlying graph is constructed
according to the [Erdős–Rényi model
](https://doi.org/10.5486%2FPMD.1959.6.3-4.12). The keyword argument `verts`
specifies the possible values for the parameter \$n\$ (i.e., how many vertices
to generate), which is generated using `rand`.
"""
function rand_matroid_er59(m; verts=nothing, rng=default_rng())
    if isnothing(verts)
        # Find the minimum vertices in a graph with `m` edges
        min_verts = ceil(Int, sqrt(2 * m) + (1 / 2))
        max_verts = ceil(Int, sqrt(2) * m)
        verts = min_verts:max_verts
    end

    n = rand(rng, verts)
    G = erdos_renyi(n, m, rng=rng)
    return GraphicMatroid(G)
end

rand_matroid_er59(V::Profile; kwds...) = rand_matroid_er59(ni(V); kwds...)


"""
    rand_matroid_er59_asym(n, m; kwds...)

Generate a random `Matroid` for each agent `n`.

See [`rand_matroid_er59`](@ref) for the supported keyword arguments.
"""
rand_matroid_er59_asym(n, m; kwds...) =
    [rand_matroid_er59(m; kwds...) for _ in 1:n]


# Type alias for a set of bitsets.
const Family = Set{SmallBitSet}


"""
    knuth_matroid(m, X)
    knuth_matroid(V::Profile, X)

Knuth's matroid construction (1975). Generates a matroid in terms of its closed
sets, given by the size of the universe `m` and a list of enlargements `X`.
"""
function knuth_matroid(m, X)
    m <= 64 || throw(DomainError(m, "number of items > 64 is not supported"))

    r::Int = 1 # r is current rank +1 due to 1-indexing.
    F::Vector{Family} = [Family([SmallBitSet()])]
    E::SmallBitSet = SmallBitSet(1:m)
    rank::Dict{SmallBitSet,Int} = Dict(SmallBitSet() => 0)

    while E ∉ F[r]
        # Initialize F[r+1].
        push!(F, Family())

        _generate_covers!(F, rank, r, E)

        # Perform enlargements.
        if r <= length(X)
            if X[r] !== nothing
                for x in X[r]
                    _add_set!(F, rank, x, r)
                end
            end
        end

        r += 1
    end

    return ClosedSetsMatroid(m, r - 1, F)
end

knuth_matroid(V::Profile, X) = knuth_matroid(ni(V), X)


"""
    rand_matroid_knu75(m; r=2:m÷2, track_indep=false, rng=default_rng())
    rand_matroid_knu75(V::Profile; ...)

Generate a random `Matroid` based on the process described in [Knuth's 1975
paper](https://doi.org/10.1016/0012-365X(75)90075-8). The matroid will have a
ground set of `1:m`.

The keyword argument `r` specifies the possible values for the target rank of
the matroid, and is by default `2:m÷2`. `track_indep` controls whether the
matroid generation should keep track of independent sets under construction.

!!! warning

    This implementation uses [`SmallBitSet`](@ref), which is backed by a
    `UInt64`, and therefore it does not support number of items `m` larger
    than 64.

"""
function rand_matroid_knu75(m; r=2:m÷2, track_indep=false, rng=default_rng())
    m <= 64 || throw(DomainError(m, "number of items > 64 is not supported"))

    rr = rand(rng, r)
    P = rand_matroid_coarsening(m, rr, rng=rng)

    if track_indep
        return _rand_matroid_knu75_full(m, P, rng=rng)
    else
        return _rand_matroid_knu75_closed(m, P, rng=rng)
    end
end

rand_matroid_knu75(V::Profile; kwds...) = rand_matroid_knu75(ni(V); kwds...)


"""
    rand_matroid_knu75_asym(n, m; kwds...)

Generate a random `Matroid` for each agent `n`.

See [`rand_matroid_knu75`](@ref) for the supported keyword arguments.
"""
rand_matroid_knu75_asym(n, m; kwds...) =
    [rand_matroid_knu75(m; kwds...) for _ in 1:n]


# Generate a random coarsening `P` that achieves the given rank `r` for a
# matroid of `m` items. The list `P = [p₀, p₁, …]` will be generated such that
# `p₀ = 0`, and `p₁, p₂, …` will be non-increasing.
function rand_matroid_coarsening(m, r; rng=default_rng())
    P = zeros(Int, r + 1)
    diff = m - r
    if diff <= 0 || r == 0
        return P
    end
    remainder = diff - sum(P)

    for i in 2:r+1
        x = rand(rng, 1:remainder)
        P[i] = x
        remainder -= x
        if remainder <= 0
            break
        end
    end

    if remainder > 0
        i = rand(rng, 2:max(r+1, 2))
        P[i] += remainder
    end

    @views sort!(P[2:end], rev=true)
    return P
end


function _rand_matroid_knu75_closed(m, P; rng=default_rng())
    r::Int = 1
    F::Vector{Family} = [Family([SmallBitSet()])]
    E::SmallBitSet = SmallBitSet(1:m)
    rank::Dict{SmallBitSet,Int} = Dict(SmallBitSet() => 0)

    while E ∉ F[r]
        r > m && @warn "Rank is larger than universe!" m r

        # Initialize F[r+1].
        push!(F, Family())

        _generate_covers!(F, rank, r, E)

        # Perform coarsening.
        if r <= length(P)
            _coarsen!(F, rank, P[r], r, E, rng=rng)
        end

        r += 1
    end

    return ClosedSetsMatroid(m, r - 1, F)
end


function _rand_matroid_knu75_full(m, P; rng=default_rng())
    r::Int = 1
    E::SmallBitSet = SmallBitSet(1:m)
    rank::Dict{SmallBitSet,Int} = Dict(SmallBitSet() => 0)

    F::Vector{Family} = [Family([SmallBitSet()])]
    I::Vector{Family} = [Family([SmallBitSet()])]

    while E ∉ F[r]
        r > m && @warn "Rank is larger than universe!" m r

        # Initialize F[r+1] and I[r+1].
        push!(F, Family())
        push!(I, Family())

        _generate_covers!(F, rank, r, E, I)

        # Perform coarsening.
        if r <= length(P)
            _coarsen!(F, rank, P[r], r, E, I, rng=rng)
        end

        r += 1
    end

    return FullMatroid(m, r - 1, F, I, nothing)
end


# Generates minimal closed sets for rank r+1 and inserts them into F[r+1],
# using the supplied insert_fn. This function should take one argument, the
# newly added set.
function _generate_covers!(
    F::Vector{Family},
    rank::Dict{SmallBitSet,Int},
    r::Int,
    E::SmallBitSet,
    I::Union{Vector{Family},Nothing}=nothing
)
    for y in F[r]
        t = setdiff(E, y)
        # Find all sets in F[r+1] that already contain y and remove excess
        # elements from t.
        for x in F[r+1]
            if issubset(y, x)
                setdiff!(t, x)
            end
            if isempty(t)
                break
            end
        end
        # Insert y ∪ a for all a ∈ t.
        while !isempty(t)
            a = minimum(t)
            _add_set!(F, rank, union(y, a), r, I)
            delete!(t, a)
        end
    end
end


# Apply the specified number of coarsenings to the matroid.
function _coarsen!(
    F::Vector{Family},
    rank::Dict{SmallBitSet,Int},
    count::Int,
    r::Int,
    E::SmallBitSet,
    I::Union{Vector{Family},Nothing}=nothing;
    rng=default_rng()
)
    for _ in 1:count
        if E ∈ F[r+1]
            return
        end
        A = rand(rng, F[r+1])
        t = setdiff(E, A)
        a = rand(rng, t)
        delete!(F[r+1], A)

        _add_set!(F, rank, union(A, a), r, I)
    end
end


function _add_set!(
    F::Vector{Family},
    rank::Dict{SmallBitSet,Int},
    x::SmallBitSet,
    r::Int,
    I::Union{Vector{Family},Nothing}=nothing
)
    done = false

    while !done
        done = true

        for y in F[r+1]
            """
            When we are comparing a set X which we are investigating whether is a
            closed set, with a set Y which we know (based on the sets that have been
            added thus far) is a closed set, we can encounter three scenarios:

            1. X ∩ Y is a closed set (it is in our rank table) of rank < r
                - Move on. If this is the case for all Y ∈ F[r+1] we add X to F[r+1].

            2. X ∩ Y is a closed set (it is in our rank table) of rank == r
                - X and Y are not closed sets. Remove Y from F[r+1] and call this
                function on X ∪ Y.

            3. We have not seen X ∩ Y before (this happens when we have closed sets
                of lower rank but similar cardinality).

                a. If |X ∩ Y| < r, we know that the rank of X ∩ Y is < r. Move on.

                b. If |X ∩ Y| >= r, we need to see if X ∩ Y ⊆ Z for some Z of lower
                    rank.
                    If not, remove Y from F[r+1] and call this function on X ∪ Y.
            """

            xy = intersect(x, y)
            xy_rank = get(rank, xy, nothing)
            if !isnothing(xy_rank)
                if xy_rank < r
                    continue
                end
            else
                if length(xy) < r
                    continue
                else
                    r´ = _check_rank(F, r, xy)
                    if !isnothing(r´)
                        rank[xy] = r´
                        continue
                    end
                end
            end

            # x ∩ y has rank > r, replace with x ∪ y.
            delete!(F[r+1], y)
            union!(x, y)
            done = false
            break
        end
    end

    push!(F[r+1], x)

    if isnothing(I)
        rank[x] = r
    else
        _add_set_callback!(I, rank, x, r)
    end
end


# Given a closed set x,
# 1. simply return if rank[x] < r (we've seen this already)
# 2. add it to I if |x| = r
# 3. recursively call this func on all x' ⊂ x st |x'| = |x| - 1
function _add_set_callback!(I::Vector{Family}, rank::Dict{SmallBitSet,Int}, x::SmallBitSet, r::Int)
    c = length(x)

    if haskey(rank, x) && rank[x] <= r
        return
    end
    if c == r
        push!(I[r+1], x)
    end
    rank[x] = r
    t = x
    while !isempty(t)
        v = setdiff(t, minimum(t))
        _add_set_callback!(I, rank, setdiff(x, minimum(union(t, v))), r)
        t = v
    end
end


function _check_rank(F::Vector{Family}, r::Int, v::SmallBitSet)
    for (i, Fi) in enumerate(@view F[1:r])
        for z ∈ Fi
            if issubset(v, z)
                return i - 1
            end
        end
    end
    return nothing
end


"""
    knuth_matroid_erect(m, enlargements)

An improved version of KMC we are also finding the independent sets and circuits
of the matroid during generation.

This version assigns the Hamming weight of all subsets of E upfront. This is
infeasible for values of n much larger than 16.
"""
function knuth_matroid_erect(m, enlargements)
    m <= 64 || throw(DomainError(m, "number of items > 64 is not supported"))

    # Initialize.
    r::Int = 1
    mask::SmallBitSet = SmallBitSet(1:m)
    mask_bits::Int = 2^m - 1
    rank::Dict{SmallBitSet,Int} = Dict()

    # Populate rank table with 100+cardinality for all subsets of E.
    rank[SmallBitSet()] = 100
    k = 1
    while (k <= mask_bits)
        for i in 0:k-1
            rank[bits_to_set(k + i)] = rank[bits_to_set(i)] + 1
        end
        k = k + k
    end

    F = [Family([SmallBitSet()])] # F[r] is the family of closed sets of rank r-1.
    I = [Family([SmallBitSet()])] # I[r] is the family of independent sets of rank r-1.
    rank[SmallBitSet()] = 0

    while mask ∉ F[r]
        push!(F, Family())
        push!(I, Family())

        # Generate minimal closed sets for rank r+1.
        for y in F[r] # y is a closed set of rank r.
            t = setdiff(mask, y) # The set of elements not in y.
            # Find all sets in F[r+1] that already contain y and remove excess elements from t.
            for x in F[r+1]
                if issubset(y, x)
                    setdiff!(t, x)
                end
            end
            # Insert y ∪ a for all a ∈ t.
            while !isempty(t)
                x = union(y, minimum(t))
                insert_set!(x, F, r, rank)
                setdiff!(t, x)
            end
        end

        # Enlarge (if any).
        if r <= length(enlargements) && enlargements[r] !== nothing
            for set in enlargements[r]
                insert_set!(set, F, r, rank)
            end
        end

        # Assign rank to sets and add independent ones to I.
        for M in F[r+1]
            mark!(M, I, r, rank)
        end

        # Next rank.
        r += 1
    end

    C = Family() # C is the set of circuits (minimal dependent sets) for M.
    k = 1
    while k <= mask_bits
        for i in 0:k-1
            ki_set = bits_to_set(k + i)
            i_set = bits_to_set(i)
            if rank[ki_set] == rank[i_set]
                push!(C, ki_set)
                unmark!(ki_set, rank[i_set] + 101, rank, mask)
            end
        end
        k += k
    end

    return FullMatroid(m, r - 1, F, I, C)
end

knuth_matroid_erect(V::Profile, enlargements) = knuth_matroid_erect(ni(V), enlargements)


# Inserts set x into F[r+1], but augments x if it is necessary to ensure no two
# sets in F[r+1] have an intersection of rank greater than r.
function insert_set!(x, F, r, rank)
    for y in F[r+1] # +1 since Julia is 1-indexed.
        if rank[intersect(x, y)] < r
            continue
        end

        # x ∩ y has rank > r, replace x and y with x ∪ y.
        delete!(F[r+1], y)
        insert_set!(union(x, y), F, r, rank)
        return
    end

    push!(F[r+1], x)
end


# Given a closed set m, sets rank[m']=r for all subsets m' ⊆ m whose rank is not
# already ≤ r, and adds m' to I if it is independent (that is, if its rank equals
# its cardinality).
function mark!(m, I, r, rank)
    if haskey(rank, m) && rank[m] <= r
        return
    end
    if rank[m] == 100 + r
        push!(I[r+1], m)
    end
    rank[m] = r
    t = m
    while !isempty(t)
        v = setdiff(t, minimum(t))
        mark!(setdiff(m, minimum(union(t, v))), I, r, rank)
        t = v
    end
end


function unmark!(m, card, rank, mask)
    if rank[m] < 100
        rank[m] = card
        t = setdiff(mask, m)
        while !isempty(t)
            v = setdiff(t, minimum(t))
            unmark!(union(m, minimum(union(t, v))), card + 1, rank, mask)
            t = v
        end
    end
end
