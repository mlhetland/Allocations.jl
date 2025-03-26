function m_check_ef_(V, A, value_; indep=false)
    N = agents(V)

    for i in N, j in N
        i !== j || continue

        val = indep ? length(bundle(A, i)) : value(V, i, A)

        if val < value_(V, i, bundle(A, j))
            return false
        end
    end

    return true
end


# General check for threshold-bases fairness measures (PROP, PROP1, PROPX, MMS),
# where each agent's bundle value should exceed some threshold. This threshold
# is computed with the supplied function thresh.
function check_thresh(V, A, thresh; indep=false)
    for i in agents(V)
        val = indep ? length(bundle(A, i)) : value(V, i, A)

        if val < thresh(V, i, A)
            return false
        end
    end

    return true
end


function check_alpha_ef_(V, A, value_; indep=false)
    N = agents(V)
    alphas = []

    for i in N, j in N
        i !== j || continue

        val = indep ? length(bundle(A, i)) : value(V, i, A)

        push!(alphas, val / value_(V, i, bundle(A, j)))
    end

    return minimum(alphas)
end

function check_alpha_thresh(V, A, thresh; indep=false)
    alphas = []
    for i in agents(V)
        val = indep ? length(bundle(A, i)) : value(V, i, A)

        push!(alphas, val / thresh(V, i, A))
    end

    return minimum(alphas)
end


check_prop(V::MatroidRank, A; indep=false) =
    check_thresh(V, A, prop; indep=indep)

check_prop1(V::MatroidRank, A; indep=false) =
    check_thresh(V, A, prop_1; indep=indep)

check_propx(V::MatroidRank, A; indep=false) =
    check_thresh(V, A, prop_x; indep=indep)

check_propx0(V::MatroidRank, A; indep=false) =
    check_thresh(V, A, prop_x0; indep=indep)


check_ef(V::MatroidRank, A; indep=false) =
    m_check_ef_(V, A, value; indep=indep)

check_ef1(V::MatroidRank, A; indep=false) =
    m_check_ef_(V, A, value_1; indep=indep)

check_efx(V::MatroidRank, A; indep=false) =
    m_check_ef_(V, A, value_x; indep=indep)

check_efx0(V::MatroidRank, A; indep=false) =
    m_check_ef_(V, A, value_x0; indep=indep)


function check_mms(V, A; indep=false)
    mmss = [mms_i(V, i) for i in agents(V)]
    for i in agents(V)
        val = indep ? length(bundle(A, i)) : value(V, i, A)

        if val < mmss[i]
            return false
        end
    end

    return true
end


function check_alpha_mms(V, A; indep=false)
    mmss = [mms_i(V, i) for i in agents(V)]
    alphas = []
    for i in agents(V)
        val = indep ? length(bundle(A, i)) : value(V, i, A)

        push!(alphas, val / mmss[i])
    end

    return minimum(alphas)
end


## ENVY-FREENESS ####################################################


function value_1(V::MatroidRank, i, S)
    if is_indep(V.Ms[i], S)
        return max(length(S) - 1, 0)
    end

    return minimum(value(V, i, setdiff(S, g)) for g in items(V))
end


# EFX and EF1 is equivalent (the least positively valued item has the same 
# value as the highest valued item).
value_x(V::MatroidRank, i, A) = value_1(V, i, A)

# If S is independent for M_i, the least valued good in S has value 1, else 0.
value_x0(V::MatroidRank, i, S) =
    is_indep(V.Ms[i], S) ? max(length(S) - 1, 0) : value(V, i, S)


## THRESHOLDS #######################################################


"""
    prop(V::MatroidRank, i, S)

The proportional share of agent `i`, that is, the value `i` places on the set 
of all goods, divided by the numder of agents.
"""
prop(V::MatroidRank, i, _) = rank(V.Ms[i]) / na(V)


"""
    prop_1(V::Profile, i, S)

The proportional share of agent `i`, minus a highest-valued item not allocated
to `i`.
"""
prop_1(V::MatroidRank, i, A) = max(prop(V, i, A) - 1, 0)


"""
    prop_x(V::Profile, i, S)

The proportional share of agent `i`, minus a least positively-valued item not allocated to `i`.
"""
prop_x(V::MatroidRank, i, A) = prop_1(V, i, A)


"""
    prop_x0(V::Profile, i, S)

The proportional share of agent `i`, minus a least-valued (possibly 0-valued)
item not allocated to `i`.
"""
prop_x0(V::MatroidRank, i, A) =
    is_closed(V.Ms[i], bundle(A, i)) ? prop_1(V, i, A) : prop(V, i, A)


"""
    function mms_i(V, i)

Finds the maximin share of agent i in the instance V.
"""
function mms_i(V::MatroidRank, i)
    M_i = V.Ms[i]
    n = na(V)

    # An initial partition into independent subsets (subjectively so for i).
    (A, _) = matroid_partition_knuth73([M_i for _ in 1:n])

    # Setup matrix D st D[j,k] v_i(A_j) - v_i(A_k) ∀ j,k ∈ [n].
    D = zeros(Integer, n, n)
    for j in 1:n
        for k in 1:n
            # v_i(A_p) = |A_p| since all sets in A are independent wrt M_i.
            D[j, k] = length(A[j]) - length(A[k])
        end
    end

    j, k = argmax(D) |> Tuple
    while D[j, k] > 1

        # By the augmentation property, ∃g ∈ A_j st A_k + g ∈ I_i.
        g = nothing
        for g´ ∈ A[j]
            if is_indep(M_i, A[k] ∪ g´)
                g = g´
                break
            end
        end

        if g === nothing
            display(D)
            display(A)
            display(V.Ms[i])
            throw("Did not find g ∈ A_$j st A_$k + g ∈ I_$i !")
        end

        # Update A.
        setdiff!(A[j], g)
        union!(A[k], g)

        # Update D.
        for l in 1:n
            D[j, l] -= 1
            D[l, j] += 1 # A_j is one smaller.
            D[k, l] += 1
            D[l, k] -= 1 # A_k is one larger.
        end

        j, k = argmax(D) |> Tuple
    end

    return minimum(length, A)
end
