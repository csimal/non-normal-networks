"""
    laplacian_coupling!(du,dv,u,v,p,t)

The laplacian coupling between two nodes `u` and `v` along the directed edge `(v,u)`.

This yields the standard diffusion on a graph (i.e. with the laplacian matrix).
"""
function laplacian_coupling!(du,dv,u,v,p,t)
    @. du = v # what goes from v to u
    @. dv = -v # comes out of v
    nothing
end

"""
    sine_coupling(du,dv,u,v,p,t)

The sine coupling between nodes `u` and `v` along the directed edge `(v,u)`.

This is the coupling used in the Kuramoto model.
"""
function sine_coupling!(du,dv,u,v,p,t)
    @. du = sin(v-u)
    @. dv = 0.0 # This has to be zero, otherwise it wouldn't work for undirected networks
    nothing
end

"""
    reaction_coupling(h, dim, f!, g!, D)

Construct a general reaction-coupling system on graph `h` with local dynamics `f!`, coupling function `g!` and coupling coefficients `D`.

The mathematical form of the resulting system is
```math
\dot{x_i} = f(x_i) 
+ D \sum_{(j,i)\in E} w_{ij} g_{in}(x_i,x_j) 
+ D \sum_{(i,j)\in E} w_{ji} g_{out}(x_i,x_j)
```
where ``w_{ij}`` is the weight of the edge from ``j`` to ``i``, and ``g_{in}``, ``g_{out}`` are respectively the coupling functions for the destination and the source of an edge.

Arguments:
* `h` the graph on which the system is defined. If it is a `SimpleWeightedGraph`, the weights are assumed to represent coupling coefficients.
* `dim` the dimension of the state space for each node. The resulting system uses a `(dim,nv(g))` to store its state, with column `i` storing the state of node `i`.
* `f!` must be an in-place function that computes the local dynamics for a single node.
* `g! = laplacian_coupling` must be of the form `g!(du,dv,u,v,p,t)` and modify `du` and `dv` with the coupling terms for the directed edge `(v,u)`. For undirected networks, this function gets called for both `(u,v)` and `(v,u)` so be careful when defining your own `g!`. The default is the laplacian coupling, yielding a reaction-diffusion system.
* `D = I` the coupling coefficients. can be any `(dim,dim)` matrix or just a scalar. This is typically a diagonal matrix.
"""
function reaction_coupling(
    h::AbstractGraph, 
    dim::Integer,
    f!::Function,
    g!::Function = laplacian_coupling!,
    D = I
    )
    di = Vector{Float64}(undef, dim)
    dj = Vector{Float64}(undef, dim)
    # Notice that we only use the adjacency matrix so this is agnostic on whether the graph is directed/weighted
    ijw = zip(findnz(adjacency_matrix(h, dir = :in))...)
    return function(du, u, p, t)
        # first compute the reaction terms
        for i in 1:nv(g)
            f!(@view du[:,i], @view u[:,i], p, t) 
        end
        for (i,j,w) in ijw
            g!(di, dj, @view u[:,i], @view u[:,j])
            mul!(@view du[i], D, di, w, 1.0) # du[i] += w * D * di
            mul!(@view du[j], D, dj, w, 1.0)
        end
    end
end