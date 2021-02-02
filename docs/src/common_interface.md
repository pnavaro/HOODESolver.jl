# SplitODEProblem

The package [DiffentialEquations.jl](https://diffeq.sciml.ai/stable/types/split_ode_types/) 
offers a interface for problems similar to the ones we are trying to solve.

To solve the Henon-Heiles system with it, you do as follows:

```
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

prob = SplitODEProblem(f1, f2, u0, tspan);
sol = solve(prob1, CNAB2(), saveat=(0:0.001:0.1));

plot(sol, vars=[3])
```

