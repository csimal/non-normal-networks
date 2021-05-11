
function _henrici(frob, λ, normalize)
    if normalize
        return sqrt( frob^2 - sum(abs.(λ).^2) ) / frob
    else
        return sqrt( frob^2 - sum(abs.(λ).^2) )
    end
end

"""
    henrici(A; normalize=true)

Compute Henrici's departure from normality for matrix A.

Henrici's DFN is defined as
```math
\sqrt{\|A\|_F - \sum_{i=1}^n |\lambda_i|^2},
```
where ``\|~\|_F`` is the Frobenius norm and ``\lambda_i`` are the eigenvalues of A.

keyword arguments:
* `normalize=true`: Whether or not to normalize by the Frobenius norm of A.
"""
function henrici(A::AbstractMatrix; normalize=true)
    λ = eigvals(A)
    frob = norm(A) # by default, calling norm on a matrix returns the Frobenius norm
    return _henrici(frob, λ, normalize)
end

henrici(A::AbstractSparseMatrix; normalize=true) = henrici(Array(A), normalize=normalize)

function henrici(g::AbstractGraph; normalize=true)
    λ = adjacency_spectrum(g, dir=:in)
    frob = norm(adjacency_matrix(g, dir=:in))
    return _henrici(frob, λ, normalize)
end