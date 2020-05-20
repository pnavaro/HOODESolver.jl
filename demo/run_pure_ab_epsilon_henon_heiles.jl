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
using Plots.PlotMeasures
function fctmain(n_tau, prec)
    setprecision(prec)
    u0=BigFloat.([90, -44, 83, 13]//100)
    tab_eps = zeros(BigFloat,8)
    epsilon=big"0.19"
    for i=1:size(tab_eps,1)
        tab_eps[i] = epsilon
        epsilon /= 10
    end
    nbmaxtest=12
    order=6
    t_max = big"1.0"
    y = ones(Float64, nbmaxtest, size(tab_eps,1) )
    x=zeros(Float64,nbmaxtest)
    ind=1
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for epsilon in tab_eps
        prob = HiOscODEProblem(henon_heiles, u0, (big"0.0",t_max), missing, A, epsilon)
        tabsol = Array{Array{BigFloat,1},1}(undef,0)
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
            println("epsilon=$epsilon\nresult=$y\nlog(result)=$(log2.(y))")
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
        Plots.savefig(p,"out/r_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_epsilon_hh_2.pdf")
        ind+= 1
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctmain(32)
