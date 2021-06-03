# Bipartite matching, using the Ford–Fulkerson algorithm (with DFS). (Algorithms
# with better asymptotic running time exist, but they may not be more efficient
# in practice, and they tend to be quite a bit more complex.) `X` should be a
# matrix where `X[i, g] ≠ 0` indicates an edge between `i` and `j`. I.e., this
# may be a boolean matrix, or any other numeric matrix where a matching along
# nonzero values (e.g., using the valuation matrix directly in the case of MNW)
# is desired. The function returns an iterator over pairs `(i => g)`,
# indicating that `i` is matched with `g`.
function bipartite_matching(X)

    n, m = size(X)

    mate, seen = zeros(Int, m), falses(m)

    function match!(i)
        for g = 1:m
            (iszero(X[i, g]) || seen[g]) && continue
            seen[g] = true
            if mate[g] == 0 || match!(mate[g])
                mate[g] = i
                return true
            end
        end
        return false
    end

    for i = 1:n
        seen .= false
        match!(i)
    end

    return ((i => g) for (g, i) in pairs(mate) if i ≠ 0)

end


bipartite_matching(X::Profile) = bipartite_matching(matrix(X))
