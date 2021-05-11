
function random_walk_laplacian(g::AbstractGraph)
    A = adjacency_matrix(g, dir=:in)
    k = reshape(sum(A, dims=2), nv(g)) # this gives the weighted degree, instead of indegree(g) which counts ingoing links
    d = [x==0.0 ? 0.0 : 1.0/x for x in k]
    return I - Diagonal(d) * A
end

"""
    laplacian(g)

Return the laplacian matrix of graph `g`.

This differs from `laplacian_matrix` in that it handles directed networks correctly for solving the heat equation.
"""
function laplacian(g::AbstractGraph)
    A = adjacency_matrix(g, dir=:in)
    d = reshape(sum(A, dims=1), nv(g)) # out-strength
    return Diagonal(d) - A
end