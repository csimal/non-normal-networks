
"""
    add_backward_loops!(g, w)

Add a backward link to each directed link in `g` with weight `w`.
"""
function add_backward_links!(g::MetaDiGraph, w::Real)
    for e in collect(edges(g))
        (u, v) = src(e), dst(e)
        if !has_edge(g, v, u)
            add_edge!(g, v, u)
            set_prop!(g, v, u, :weight, w)
        end
    end
    return g
end