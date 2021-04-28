
"""
    add_reciprocal_links!(g, p, w)

Add reciprocal links to a weighted digraph, with probability `p` and weights sampled from `w`.

If `w` is a number, then it is assigned as the weight for every created link.
"""
function add_reciprocal_links!(g::MetaDiGraph, p::Real, w::Distribution = Uniform(0.0,1.0))
    for e in edges(g)
        (u, v) = src(e), dst(e)
        if !has_edge(g, v, u) && rand() < p
            add_edge!(g, v, u)
            set_prop!(g, v, u, :weight, rand(w))
        end
    end
    return g
end

function add_reciprocal_links!(g::MetaDiGraph, p::Real, w::Real)
    for e in edges(g)
        (u,v) = src(e), dst(e)
        if !has_edge(g, v, u) && rand() < p
            add_edge!(g, v, u)
            set_prop!(g, v, u, :weight, w)
        end
    end
    return g
end

"""
    add_reciprocal_links(g, p, w)

Add reciprocal links to a weighted digraph, with probability `p` and weights sampled from `w`.

This function creates a copy of `g`. See [add_reciprocal_links!](@ref) for the in-place version.
"""
add_reciprocal_links(g, p, w=Uniform(0.0,1.0)) = add_reciprocal_links!(copy(g), p, w)

add_reciprocal_links!(p::Real, w=Uniform(0.0,1.0)) = g -> add_reciprocal_links!(g, p, w)
add_reciprocal_links(p::Real, w=Uniform(0.0,1.0)) = g -> add_reciprocal_links(g, p, w)

