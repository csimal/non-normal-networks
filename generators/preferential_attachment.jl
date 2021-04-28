
"""
    preferential_attachment!(g, n, k)

Add `n` vertices to graph `g` and link them existing vertices by preferential attachment.

`k` can be
* a positive integer, in which case up to `k` links are created
* a discrete distribution, in which case the number of links to create is sampled from that distribution (truncated between 0 and the current number of vertices)
"""
function preferential_attachment!(g, n::Integer, k::Integer)
    N = nv(g)
    for i in 1:n
        # sample up to k unique neighbors proportionally to the in-degree
        targets = sample(vertices(g), Weights(indegree(g)), min(k, nv(g)), replace=false)
        add_vertex!(g)
        for j in targets
            add_edge!(g, N+i, j)
        end
    end
    return g
end

function preferential_attachment!(g, n::Integer, k::DiscreteDistribution)
    N = nv(g)
    for i in 1:n
        t = truncated(k, 1, nv(g))
        targets = sample(vertices(g), Weights(indegree(g)), rand(t), replace=false)
        add_vertex!(g)
        for j in targets
            add_edge!(g, N+i, j)
        end
    end
    return g
end

"""
    preferential_attachment(g, n, k)

Add `n` vertices to graph `g` and link them existing vertices by preferential attachment.

`k` can be
* a positive integer, in which case up to `k` links are created
* a discrete distribution, in which case the number of links to create is sampled from that distribution (truncated between 0 and the current number of vertices)

This function creates a copy of `g`. For the in-place version, see [preferential_attachment!](@ref)
"""
preferential_attachment(g, n, k) = preferential_attachment!(copy(g), n, k)

preferential_attachment!(n::Integer, k) = g -> preferential_attachment!(g, n, k)
preferential_attachment(n::Integer, k) = g -> preferential_attachment(g, n, k)
