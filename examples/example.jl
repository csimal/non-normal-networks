using DifferentialEquations
using Plots
using MATLAB
using LightGraphs
using LinearAlgebra
using LaTeXStrings

g = watts_strogatz(10, 4, 0.1)
L =  - laplacian_matrix(g) # we negate L now for the heat equation

f! = function(du, u, p, t)
    mul!(du, L, u) # This is equivalent to du = L*u, but it doesn't allocate memory for the result of L*u, instead it puts the result directly into du
end

u0 = randn(10) # we set the initial condition to be normally distributed
tspan = (0.0, 2.0)

# This constructs your ODE
prob = ODEProblem(f!, u0, tspan)

# and this solves it. Tsit5() is a 4-5 order Runge-Kutta method with adaptive timestep
sol = solve(prob,Tsit5())

plot(sol) # You can just plot the solution

f(t, u...) = (t, norm(u))
# You can also plot a function of the solution. The 0 index here is for time
plot(sol, vars=(f, 0:10...), label=L"\| u \|")

# Let's use Matlab to plot this instead
# expressions enclosed in $() are converted to Matlab format
mat"""
figure
plot($(sol.t), $(map(norm, sol.u)))
legend({'\$\\| u \\| \$'}, 'Interpreter', 'latex')
"""