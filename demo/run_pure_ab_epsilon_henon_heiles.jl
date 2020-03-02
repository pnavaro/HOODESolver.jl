#=
testODE:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-13
=#
include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using LinearAlgebra
using Plots
using Random
function fctmain(n_tau)
    seed=123456
    Random.seed!(seed)
    u0=rand(BigFloat,4)
    println("seed = $seed")
    tab_eps = zeros(BigFloat,14)
    epsilon=big"0.8"
    for i=1:14
        tab_eps[i] = epsilon
        epsilon /= 2^i
    end
    nbmaxtest=8
    order=6
    ordprep=8
    t_max = big"1.0"
    y = ones(Float64, nbmaxtest, size(tab_eps,1) )
    x=zeros(Float64,nbmaxtest)
    ind=1
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for epsilon in tab_eps
        parphi = PreparePhi(epsilon, n_tau, A, henon_heiles)
        println("prepareU0 eps=$epsilon n_tau=$n_tau")
        @time par_u0 = PrepareU0(parphi, ordprep, u0)
        nb=100*2^nbmaxtest
        pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
        solref =  twoscales_pure_ab(pargen, only_end=true)
        eps_v = convert(Float32,epsilon)
        println("epsilon = $eps_v solref=$solref")
        nb = 100
        indc =1
        labels=Array{String,2}(undef, 1, ind)  
        while indc <= nbmaxtest
            pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
            sol =  twoscales_pure_ab(pargen, only_end=true)
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            diff=solref-sol
            x[indc] = 1.0/nb
            println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
            y[indc,ind] = norm(diff,Inf)
            println("epsilon=$epsilon result=$y")
            nb *= 2
            indc += 1
        end
        for i=1:ind
            labels[1,i] = " epsilon,order=$(convert(Float32,tab_eps[i])),$order "
        end
        gr()
        p=Plots.plot(
    x,
    view(y,:,1:ind),
    xlabel="delta t",
    xaxis=:log,
    ylabel="error",
    yaxis=:log,
    legend=:topleft,
    label=labels,
    marker=2
)
        
        prec_v = precision(BigFloat)
        Plots.savefig(p,"out/r_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_epsilon_hh_1.pdf")
        ind+= 1
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctmain(32)
