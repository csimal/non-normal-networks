
"""
    trophic_level(g)

Compute the trophic level of each node in a directed graph.
"""
function trophic_level(g::AbstractGraph)
    if all(x -> x!=0, indegree(g))
        @warn "No nodes with zero indegree"
        return ones(nv(g))
    end
    L = random_walk_laplacian(g)
    return L\ones(nv(g))
end

function trophic_difference(g::AbstractGraph, h)
    return [h[dst(e)]-h[src(e)] for e in edges(g)]
end

"""
    trophic_incoherence(g)

Return the trophic incoherence of a directed graph.

A value of zero indicates a maximally coherent digraph.
"""
function trophic_incoherence(A::AbstractGraph)
    h = trophic_level(g)
    x = trophic_difference(g, h)
    return std(x, mean=1.0)
end