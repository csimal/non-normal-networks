using Plots
using LinearAlgebra

f(t, u...) = (t, norm(u))
f1(t, u...) = (t, sum(u))

g1 = path_graph(10)
g2 = path_digraph(10)

t = 10.0
x₀ = zeros(10)
x₀[1] = 1.0

prob1 = ODEProblem(nd_heat(g1), x₀, (0.0,t))
sol1 = solve(prob1, Tsit5())
plot(sol1)

prob2 = ODEProblem(laplacian_operator(g1), x₀, (0.0,t))
sol2 = solve(prob2, Tsit5())
plot(sol2)

prob3 = ODEProblem(nd_heat(g2), x₀, (0.0,t))
sol3 = solve(prob3, Tsit5())
plot(sol3)

prob4 = ODEProblem(laplacian_operator(g2), x₀, (0.0,t))
sol4 = solve(prob4, Tsit5())
plot(sol4)

prob5 = ODEProblem(nd_heat(reverse(g2)), x₀, (0.0,t))
sol5 = solve(prob5, Tsit5())
plot(sol5)