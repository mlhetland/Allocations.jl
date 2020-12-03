# Bipartite mathcing using an LP solver, for simplicity (as we're using a MIP
# solver elsewhere, anyway). `X` should be a boolean matrix, where `X[i, j]`
# indicates an edge between `i` and `j`. The `solver` argument should be JuMP
# optimizer factory.
function bipartite_matching(X, solver=MIP_SOLVER)

    model = Model(solver)

    n, m = size(X)

    @variable(model, x[1:n, 1:m] >= 0)

    for i = 1:n, j = 1:m
        X[i, j] || fix(x[i, j], 0, force=true)
    end

    for i = 1:n
        @constraint(model, sum(x[i, j] for j = 1:m) <= 1)
    end

    for j = 1:m
        @constraint(model, sum(x[i, j] for i = 1:n) <= 1)
    end

    @objective(model, Max, sum(x))

    optimize!(model)

    @assert termination_status(model) == MOI.OPTIMAL

    return JuMP.value.(x) .≈ 1.0

end
