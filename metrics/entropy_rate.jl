

function entropy_rate(T::AbstractMatrix; tol=1e-06, nmax=1000)
    N = size(T)[1]
    q₀ = ones(Float64, N) ./ N
    q₁ = T * q₀
    n = 0
    while norm(q₀-q₁) > tol && n < nmax
        q₀ .= q₁
        mul!(q₁, T, q₀)
        n += 1
    end
    h = 0.0
    for j in 1:N, i in 1:N
        h -= T[i,j] * q₁[j] * log(T[i,j])
    end
    return h
end