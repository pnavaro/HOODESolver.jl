#=
testODE:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-13
=#
include("../src/interface.jl")
using LinearAlgebra
using Plots
using Plots.PlotMeasures
function fctMain(n_tau)

 #   u0 =[big"0.12345678", big"0.1209182736", big"0.1290582671", big"0.1239681094" ]
    tab_eps = zeros(BigFloat,7)
    epsilon=big"0.15"
    for i=1:size(tab_eps,1)
        tab_eps[i] = epsilon
        epsilon /= 10
    end
    u0 = BigFloat.([-34//100, 78//100, 67//100, -56//10])
    B = BigFloat.([12//100 -78//100 91//100 34//100
        -45//100 56//100 3//100 54//100
        -67//100 09//100 18//100 89//100
        -91//100 -56//100 11//100 -56//100])
    alpha =  BigFloat.([12//100, -98//100, 45//100, 26//100])
    beta =  BigFloat.([-4//100, 48//100, 23//100, -87//100])
    fct = (u,p,t) -> B*u + t*p[1] +p[2]
    t_max = big"1.0"
    nbMaxTest=12
    order=6
    y = ones(Float64, nbMaxTest, size(tab_eps,1) )
    x=zeros(Float64,nbMaxTest)
    ind=1
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for epsilon in tab_eps
        prob = HiOscODEProblem(fct,u0, (big"0.0",t_max), (alpha, beta), A, epsilon, B)
        eps_v = convert(Float32,epsilon)
        println("epsilon = $eps_v")
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, ind)  
        while indc <= nbMaxTest
            res = solve(prob, nb_tau=n_tau, order=order, nb_t=nb, dense=false)
            sol = res[end]
            solref = getexactsol(res.par_u0.parphi, u0, t_max)
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            diff=solref-sol
            x[indc] = t_max/nb
            println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
            ni = norm(diff,Inf)
            y[indc,ind] = (ni < 1) ? ni : NaN
            println("epsilon=$epsilon result=$y")
            println("epsilon=$epsilon reslog2=$(log2.(y))")
            nb *= 2
            indc += 1
        end
        for i=1:ind
            labels[1,i] = " ε = $(convert(Float32,tab_eps[i])) "
        end
        p=Plots.plot(
    x,
    view(y,:,1:ind),
    xlabel="Δt",
    xaxis=:log,
    ylabel="error",
    yaxis=:log,
    legend=:bottomright,
    label=labels,
    marker=2
)    
pp=Plots.plot(
    x,
    view(y,:,1:ind),
    xlabel="Δt",
    xaxis=:log,
    ylabel="error",
    yaxis=:log,
    legend=:bottomright,
    label=labels,
    marker=2,
    bottom_margin=30px,
    left_margin=60px
)    
        prec_v = precision(BigFloat)
        Plots.savefig(p,"out/r10p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon.pdf")
        Plots.savefig(p,"out/r10p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon.png")
        Plots.savefig(pp,"out/r10pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon.pdf")
        Plots.savefig(pp,"out/r10pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon.png")
        ind+= 1
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctMain(32)
