# SplitODEProblem

The package [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/types/split_ode_types/) 
offers a interface for problems similar to the ones we are trying to solve.

To solve the Henon-Heiles system, you can :

```@setup 4
using Plots
```

```@example 4
using DifferentialEquations

epsilon = 0.002
A = [ 0 0 1 0 ;
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

tspan = (0.0, 0.1)

u0 = [0.55, 0.12, 0.03, 0.89]

prob1 = SplitODEProblem(f1, f2, u0, tspan);
sol1 = solve(prob1, CNAB2(), dt=0.0001);

```

With our method we need to give the value of `epsilon` so we can't use the 
`DiffEqArrayOperator` so we define a new type:

```@example 4
using HOODESolver

linop = LinearHOODEOperator(epsilon, A)
prob2 = SplitODEProblem(linop, f2, u0, tspan)
sol2 = solve(prob2, HOODEAB(), dt=0.01)

plot(sol1, vars=[3])
plot!(sol2, vars=[3])
plot!(sol2.t, getindex.(sol2.u, 3), m=:o)
```
