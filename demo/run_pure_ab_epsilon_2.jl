#=
testODE:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-13
=#
include("../src/twoscales_pure_ab.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
using Random
function twoscales_solve( par_u0::PrepareU0, order, t, nb)
end
function fctMain(n_tau)
    seed=123321
    Random.seed!(seed)
    u0 = rand(BigFloat, 4)
    B = 2rand(BigFloat, 4, 4) - ones(BigFloat, 4, 4)
    println("seed = $seed B=$B")
    nbeps=48
    nbmaxtest=3
    paseps = big"10.0"^0.125
    ordmin=3
    ordmax=14
    ordpas = 1
    t_max = big"1.0"
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for order = ordmin:ordpas:ordmax
        y = ones(Float64, nbeps, nbmaxtest)
        x=zeros(Float64,nbeps)
        for i_eps=1:nbeps
            epsilon = big"1.0"/paseps^i_eps
            x[i_eps] = epsilon
            fct = u -> B*u
            parphi = PreparePhi(epsilon, n_tau, A, fct, B)
            println("prepareU0 eps=$epsilon n_tau=$n_tau")
            @time par_u0 = PrepareU0(parphi, order+2, u0)
            solref = getexactsol(parphi, u0, t_max)
            eps_v = convert(Float32,epsilon)
            println("epsilon = $eps_v solref=$solref")
            for indc=1:nbmaxtest
                nb = 10*10^indc
                @time parGen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
                @time sol = twoscales_pure_ab(parGen, only_end=true)
                println("solref=$solref")
                println("nb=$nb sol=$sol")
                diff=solref-sol
                println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
                y[i_eps, indc] = norm(diff,Inf)
                println("epsilon=$epsilon result=$y")
            end
        end
        labels=Array{String,2}(undef, 1, nbmaxtest)  
        for i=1:nbmaxtest
            nb = 10*10^i
            labels[1,i] = " delta t,order=$(1/nb), $order "
            mindiff = minimum(y[:,i])
            bornediff = min(mindiff*1e20,1)
            for j=1:nbeps
                if y[j,i] >= bornediff
                    y[j,i] = NaN
                end
            end
        end
        gr()
        p=Plots.plot(
    x,
    y,
    xlabel="epsilon",
    xaxis=:log,
    ylabel="error",
    yaxis=:log,
    legend=:topleft,
    label=labels,
    marker=2
)
        prec_v = precision(BigFloat)
        Plots.savefig(p,"out/r_$(prec_v)_$(order)_$(n_tau)_epsilon_lgn_2.pdf")
    end
end
fctMain(32)
