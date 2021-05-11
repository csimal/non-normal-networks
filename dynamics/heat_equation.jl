
"""
    laplacian_operator(g, σ=1.0)

Return the laplacian operator of graph `g`, as a function that can be passed to `ODEProblem`.

An optional scaling factor `σ` can be provided.

Keyword arguments:
* `overwrite=true`: whether the returned function should overwrite it's input. If set to false, the function will instead add its result it.
"""
function laplacian_operator(g::AbstractGraph, σ = 1.0; overwrite=true)
    L = -σ * laplacian(g)
    if overwrite
        return function(du,u,p,t)
            mul!(du, L, u)
        end
    else
        return function(du,u,p,t)
            mul!(du, L, u, 1.0, 1.0)
        end
    end
end



# NetworkDynamics functions. These do not handle directed networks correctly
function heat_edge!(e, v_s, v_d, p, t)
    e .= v_s - v_d
    nothing
end

function heat_vertex!(dv, v, edges, p, t)
    dv .= 0.0
    for e in edges
        dv .+= e
    end
    nothing
end

nd_heat_edge = StaticEdge(f! = heat_edge!, dim = 1)
nd_heat_vertex = ODEVertex(f! = heat_vertex!, dim = 1)
nd_heat(g) = network_dynamics(nd_heat_vertex, nd_heat_edge, g)
