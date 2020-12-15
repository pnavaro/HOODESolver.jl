module BenchParticle

using Random
using BenchmarkTools
using HOODESolver

function fct(u, p, t)
    s1, c1 = sincos(u[1]/2)
    s2, c2 = sincos(u[2])
    s3, c3 = sincos(u[3])
    return [0, 0, u[6], c1*s2*s3/2, s1*c2*s3, s1*s2*c3]
end

Random.seed!(561909)
u0 = 2rand(BigFloat,6)-ones(BigFloat,6)
epsilon=big"1e-7"
ordmax=9
A = [0 0 0  1 0 0; 
     0 0 0  0 1 0;
     0 0 0  0 0 0; 
     0 0 0  0 1 0; 
     0 0 0 -1 0 0; 
     0 0 0  0 0 0]

t_0 = big"0.0"
t_max = big"1.0"
prob = HOODEProblem(fct, u0, (t_0, t_max), missing, A, epsilon)
nb = 10^4
n_tau = 32
suite["solve"] = @benchmarkable solve(prob, nb_t=nb, order=ordmax, 
                                      getprecision=false, nb_tau=n_tau, dense=false)

end

BenchParticle.suite
