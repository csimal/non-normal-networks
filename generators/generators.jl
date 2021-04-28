using Distributions
using LightGraphs
using SimpleWeightedGraphs
using StatsBase

include("preferential_attachment.jl")
include("reciprocal_links.jl")
include("edges_operators.jl")

"""
    add_weights(g, w)

Add weights sampled from the distribution `w` to the directed edges of `g`.
"""
function add_weights(g, w::Distribution = Uniform(0.0,1.0))
    h = MetaDiGraph(g)
    for e in edges(h)
        set_prop!(h, e, :weight, rand(w))
    end
    return h
end

