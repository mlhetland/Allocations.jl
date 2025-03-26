"""
    yankee_swap(V::MatroidRank)

Viswanathan and Zick's Yankee Swap algorithm (2022) for matroid-rank valuations. The agents are prioritized according to index.
"""
function alloc_yankee_swap_vz22(V::MatroidRank)
    n = na(V)
    m = ni(V)

    # Randomly prioritize the agents.
    agents = shuffle(1:n)

    # Agent "0" (n+1) has a corresponding zero matroid.
    Ms_ = [V.Ms..., ZeroMatroid(m)]

    A = Allocation(n + 1, m)
    give!(A, n + 1, 1:m) # The bundle of unallocated items.
    flag = falses(n)

    D = exchange_graph(Ms_, A)

    while false in flag
        # The agents whose bundle can still improve.
        T = [i for i in agents if flag[i] == false]

        # Find the agents in T with minimim value.
        T_vals = [(i, length(bundle(A, i))) for i in T]
        min_val = minimum(last, T_vals)
        T_ = [i for (i, v) in T_vals if v == min_val]

        # The highest priority agent with minimum value.
        i = T_[1]

        # The goods for which i has positive marginal value.
        F_i = [g for g in 1:m if is_indep(V.Ms[i], bundle(A, i) ∪ g)]

        #  a shortest path from F_i to an unallocated good.
        A_0 = [g for g in 1:m if owner(A, g) == n + 1]
        transfer_path = find_shortest_path(D, F_i, A_0)

        # Transfer if path exists.
        if transfer_path !== nothing
            transfer!(Ms_, D, A, i, transfer_path)
        else
            flag[i] = true
        end
    end

    return A
end


"""
    alloc_eit_bciz21(V::MatroidRank)

The Envy-Induced Transfer algorithm, Algorithm 1 in Benabbou, Chakraborty, Igarashi and Zick (2021) computes a MAX-USW, EF1 
allocation.
"""
function alloc_eit_bciz21(V::MatroidRank; partition=nothing)
    n = na(V)
    m = ni(V)

    if partition === nothing
        # Compute a clean, MAX-USW allocation.
        (partition, _junk) = matroid_partition_knuth73(V.Ms)
    end

    A = Allocation(n, m)
    for (i, bundle) in enumerate(partition)
        give!(A, i, bundle)
    end

    # Envy table envy[i,j] holds i's envy towards j, v_i(A_j) - v_i(A_i).
    envy = zeros(Int, n, n)
    for i in 1:n, j in 1:n
        # We use length when we know the bundles are independent.
        envy[i, j] = value(V, i, bundle(A, j)) - length(bundle(A, i))
    end

    # While there are agents i, j st i envies j more than 1...
    i, j = argmax(envy) |> Tuple
    while envy[i, j] > 1
        # Find item in A_j with marginal gain for i.
        for g in bundle(A, j)
            if marginal_value(V, A, i, g) == 1
                # Envy-induced transfer:
                deny!(A, j, g)
                give!(A, i, g)

                # Update D.
                for k in 1:n
                    envy[i, k] = value(V, i, bundle(A, k)) - length(bundle(A, i))
                    envy[k, i] = value(V, k, bundle(A, i)) - length(bundle(A, k))
                    envy[j, k] = value(V, j, bundle(A, k)) - length(bundle(A, j))
                    envy[k, j] = value(V, k, bundle(A, j)) - length(bundle(A, k))
                end

                break
            end
        end

        i, j = argmax(envy) |> Tuple
    end

    return A
end


"""
    alloc_algmms_bv21(V::MatroidRank)

Barman and Verma's ALGMMS (2021), which computes a utilitarian
social welfare-maximizing and MMS-fair allocation.
"""
function alloc_algmms_bv21(V::MatroidRank; mmss=nothing, partition=nothing, junk=nothing)
    n = na(V)
    m = ni(V)

    if partition === nothing
        # Compute a clean, MAX-USW allocation.
        (partition, junk) = matroid_partition_knuth73(V.Ms)
    end

    A = Allocation(n, m)
    for (i, bundle) in enumerate(partition)
        give!(A, i, bundle)
    end

    if mmss === nothing
        # Compute MMS of each agent.
        mmss = [mms_i(V, i) for i in 1:n]
    end

    S_less = Set([i for i in 1:n if value(V, i, A) < mmss[i]])
    S_more = Set([i for i in 1:n if value(V, i, A) > mmss[i]])

    D = exchange_graph(V.Ms, A)

    while length(S_less) > 0
        # i is an agent with less than their maximin share.
        i = popfirst!(collect(S_less))

        # The goods for which i has positive marginal value.
        F_i = [g for g in 1:m if is_indep(V.Ms[i], bundle(A, i) ∪ g)]
        A_more = reduce(∪, [bundle(A, j) for j in S_more])

        transfer_path = find_shortest_path(D, F_i, A_more)
        @assert transfer_path !== nothing

        j = owner(A, transfer_path[end]) # The losing agent.

        transfer!(V.Ms, D, A, i, transfer_path)

        # Only i and j have received a new value.
        for k in [i, j]
            if value(V, k, A) < mmss[k]
                push!(S_less, k)
            else
                setdiff!(S_less, k)
            end
            if value(V, k, A) > mmss[k]
                push!(S_more, k)
            else
                setdiff!(S_more, k)
            end
        end
    end

    # Give agent 1 any unallocated items (these are 0-valued by everyone).
    give!(A, 1, junk)
    return A
end
