
"""
    rewire_edges(g, p)

Rewire each edge of `g` with probability `p`.

When `g` is a cyclic k-lattice, this essentially produces Watts-Strogatz networks.
"""
function rewire_edges(g::T, p::Real) where T <: AbstractGraph
    h = typeof(g)(nv(g))
    for e in edges(g)
        (u,v) = src(e), dst(e)
        # this is somewhat inefficient, but at least it's correct
        rewirings = append!(
            [(u,w) for w in setdiff(vertices(g), outneighbors(h, u)) if w != v], 
            [(w,v) for w in setdiff(vertices(g), inneighbors(h, v)) if w != u]
            )
        if !isempty(rewirings) && rand() < p
            (u, v) = rand(rewirings)
        end
        add_edge(h, u, v)
    end
    return h
end

rewire_edges(p::Real) = g -> rewire_edges(g,p)

"""
    drop_edges!(g, p)

Remove each edge in `g` with probability `p`.
"""
function drop_edges!(g, p::Real)
    for e in collect(edges(g))
        if rand() < p
            rem_edge(g,p)
        end
    end
    return g
end

"""
    drop_edges(g, p)

Remove each edge in `g` with probability `p`.

This function creates a copy of `g`. See [drop_edges!](@ref) for the in-place version.
"""
drop_edges(g, p::Real) = drop_edges!(copy(g), p)

drop_edges!(p::Real) = g -> drop_edges!(g,p)
drop_edges(p::Real) = g -> drop_edges!(g,p)

