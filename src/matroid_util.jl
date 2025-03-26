# Constructs an exchange graph over `m` matroids (E, I_i) and `n` subsets
# A_i ⊆ E, given as an n x m BitMatrix A, where A[i,j] = 1 iff element 
# j ∈ A_i. It is assumed each A_i is independent in M_i.
# The exchange graph is a graph whose vertices are the elements of E, and
# contains the edge (i,j) iff rank_i(S_i - i + j) = rank_i(S_i), where S_i is
# the set that contains i and rank_i the rank function for the corresponding
# matroid. In other words, it contains an edge (i,j) iff i can be replaced by j
# with no loss in rank.
function exchange_graph(Ms::Vector{T}, A::Allocation; all_indep=true) where {T<:Matroid}
    m = ni(A)
    D = SimpleDiGraph(m)

    # Checking for each i,j whether element ei can be replaced with ej.
    for gi in 1:m, gj in setdiff(1:m, gi)
        if !owned(A, gi)
            continue
        end
        i = owner(A, gi)

        if all_indep
            # Check if A_i - ei + ej is independent in Mi.
            if is_indep(Ms[i], setdiff(bundle(A, i), gi) ∪ gj)
                add_edge!(D, gi, gj)
            end
        else
            if rank(Ms[i], bundle(A, i)) == rank(Ms[i], setdiff(bundle(A, i), gi) ∪ gj)
                add_edge!(D, gi, gj)
            end
        end
    end

    return D
end


# Finds a shortest path [s,...,t], where s ∈ from and t ∈ to.
function find_shortest_path(D, from, to)
    X = intersect(from, to)
    if length(X) > 0
        return [X[1]]
    end

    ds = dijkstra_shortest_paths(D, from)
    paths = []

    for g in to
        path = [ds.parents[g], g]

        while path[1] ∉ from
            if path[1] == 0
                @goto skip
            end
            pushfirst!(path, ds.parents[path[1]])
        end

        push!(paths, path)
        @label skip
    end

    if length(paths) == 0
        return nothing
    end
    return argmin(length, paths)
end


# Creates a new allocation from A, where the goods have between
# transferred along path, path[1] ending up in agent i's bundle.
# It is assumed the allocation in A is clean, ie. each A_i is
# independent in M_i.
# Returns the augmented allocation and an updated exchange graph D.
function transfer!(Ms, D, A, i, path; all_indep=true)
    # At every iteration, x receives the next good in the path.
    x = i
    for g in path
        # y is the current owner, who loses g.
        y = owner(A, g)
        deny!(A, y, g)
        give!(A, x, g)

        # Recalculate neighbors of g in D.
        for g_ in vertices(D)
            # Check if A_x - g + g_ is independent in Mi.
            if all_indep
                if is_indep(Ms[x], setdiff(bundle(A, x), g) ∪ g_)
                    add_edge!(D, g, g_)
                else
                    rem_edge!(D, g, g_)
                end

            else
                if rank(Ms[x], bundle(A, x)) == rank(Ms[x], setdiff(bundle(A, x), g) ∪ g_)
                    add_edge!(D, g, g_)
                else
                    rem_edge!(D, g, g_)
                end
            end
        end

        x = y
    end
end


# Knuth's 1973 Matroid Partitioning algorithm for partitioning a set into
# subsets independent in various given matroids. Knuth's description: Given k
# matroids Ms = [M1, ..., Mk] on the same ground set E, the algorithm finds a
# k-partitioning S = [S1, ..., Sk] of the elements of E such that Sj is
# independent in matroid Mj, and nj <= |Sj| <= nj´, for given limits nj and nj´. 
# Returns the tuple (S, X), where X is the set of any elements that didn't fit
# into any independent set. This implementation drops the upper limit nj´ for
# each element j (implicitly it is infinity for all matroids). Supply nj in
# array lims (lims[j] = nj).
function matroid_partition_knuth73(Ms, floors=nothing)
    n = Ms[1].n
    k = length(Ms)
    S0 = Set(1:n) # The unallocated items.
    S = [Set() for _ in 1:k] # The partition-to-be.
    color = Dict(x => 0 for x in 1:n) # color[x] = j iff x ∈ S[j].
    for y in 1:k
        color[-y] = y
    end # -y is the 'standard' element of color y.
    succ = [0 for _ in 1:n]

    floors = floors === nothing ? [0 for _ in 1:k] : floors

    # Ensure every part gets at least its lower limit.
    for j in 1:k, i in 1:floors[j]
        augment!(j, n, k, Ms, S, S0, succ, color)
    end

    # Allocate the rest.
    while S0 != Set()
        X = augment!(0, n, k, Ms, S, S0, succ, color)

        if length(X) != 0
            return (S, X)
        end
    end

    return (S, Set())
end


function augment!(r, n, k, Ms, S, S0, succ, color)
    for x in 1:n
        succ[x] = 0
    end

    A = Set(1:n)
    B = r > 0 ? Set(-r) : Set(-j for j in 1:k)

    while B != Set()
        C = Set()
        for y ∈ B, x ∈ A
            j = color[y]

            if x ∉ S[j] && is_indep(Ms[j], x ∪ setdiff(S[j], y))
                succ[x] = y
                A = setdiff(A, x)
                C = C ∪ x
                if color[x] == 0

                    # x is uncolored - transfer:
                    while x ∈ 1:n
                        y = succ[x]
                        j = color[x]

                        if j == 0
                            setdiff!(S0, x)
                        else
                            setdiff!(S[j], x)
                        end

                        j = color[y]
                        S[j] = S[j] ∪ x
                        color[x] = j
                        x = y
                    end

                    return Set()
                end
            end
        end
        B = C
    end

    # We did not find a transfer path to 0.
    return setdiff(A, reduce(∪, S))
end
