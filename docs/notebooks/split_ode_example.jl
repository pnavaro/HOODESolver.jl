# -*- coding: utf-8 -*-
using DiffEqOperators, LinearAlgebra, DifferentialEquations
using ApproxFun, Sundials, Plots

# +
ϵ = 0.1
A = DiffEqArrayOperator( [ 0 0 1 0 ;
      0 0 0 0 ;
     -1 0 0 0 ;
      0 0 0 0 ] .+ ϵ)

function h(du, u, p, t)
    du[1] = 0
    du[2] = u[4]
    du[3] = 2*u[1]*u[2]
    du[4] = -u[2] - u[1]^2 + u[2]^2 
end

# +
epsilon= 0.0001

tspan = (0.0, 3.0)

u0 = [0.55, 0.12, 0.03, 0.89]
# -

# Define u' = Au + f
prob = SplitODEProblem(A, h, u0, (0.0,10.0), similar(u0));

sol = solve(prob, KenCarp4())


plot(sol)



