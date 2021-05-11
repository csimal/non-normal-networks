using LightGraphs
using LinearAlgebra





function hierarchy(g)
    leaders = [i for i in vertices(g) if outdegree(g, i) == 0]
    if isempty(leaders)
        @warn "No nodes with outdegree 0"
        return zeroes(Int, nv(g))
    end
    return dijkstra_shortest_paths(reverse(g), leaders).dists
end