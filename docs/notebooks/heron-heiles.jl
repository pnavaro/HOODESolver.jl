# -*- coding: utf-8 -*-
using DiffEqOperators
using DifferentialEquations
using HOODESolver
using Plots

# +
epsilon = 0.002
A = Float64[ 0 0 1 0 ;
             0 0 0 0 ;
            -1 0 0 0 ;
             0 0 0 0 ]

f1 = DiffEqArrayOperator( A ./ epsilon)

function f2(du, u, p, t)
    du[1] = 0
    du[2] = u[4]
    du[3] = 2*u[1]*u[2]
    du[4] = -u[2] - u[1]^2 + u[2]^2 
end

# +
tspan = (0.0, 0.1)

u0 = [0.55, 0.12, 0.03, 0.89]
# -

# Define u' = Au + f
prob1 = SplitODEProblem(f1, f2, u0, tspan);

sol1 = solve(prob1, CNAB2(), dt=0.0001);


prob2 = HOODEProblem(h, u0, tspan, nothing, A, ϵ);


sol2 = solve(prob2, nb_t=10);

t = LinRange(tspan[1], tspan[2], length(sol2))
plot(t, getindex.(sol2.u,3), m=:o)
plot!(sol1, vars=[3])




