
"""
    add_backward_links!(g, w)

Add a backward link to each directed link in `g` with weight `w`.
"""
function add_backward_links!(g::SimpleWeightedDiGraph, w::Real)
    for e in collect(edges(g))
        (u, v) = src(e), dst(e)
        if !has_edge(g, v, u)
            add_edge!(g, v, u, w)
        end
    end
    return g
end

add_backward_links(g, w) = add_backward_links!(copy(g), w)

add_backward_links!(w::Real) = g -> add_backward_links!(g,w)
add_backward_links(w::Real) = g -> add_backward_links(g,w)
