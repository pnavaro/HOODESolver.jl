#=
testODE:
- Julia version: 
- Author: ymocquar
- Date: 2019-11-13
=#
include("../src/interface.jl")
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
        nb = 60*2^nbmaxtest
        res = solve(prob, nb_tau=n_tau, order=order, nb_t=nb, dense=false)
        par_u0=res.par_u0
        solRef = res[end]          
        nb = 100
        indc =1
        labels=Array{String,2}(undef, 1, ind)  
        while indc <= nbmaxtest
            res = solve(prob, nb_tau=n_tau, order=order, nb_t=nb,par_u0=par_u0, dense=false)
            sol = res[end]
            nm = norm(sol - solRef, Inf)
            y[indc,ind] = (nm < 1) ? nm : NaN  
            x[indc] = 1.0/nb
            nb *= 2
            indc += 1
        end
        for i=1:ind
            labels[1,i] = "  ε = $(convert(Float32,tab_eps[i])) "
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
        Plots.savefig(p,"out/r2p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon_hh_2.pdf")
        Plots.savefig(p,"out/r2p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon_hh_2.png")
        Plots.savefig(pp,"out/r2pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon_hh_2.pdf")
        Plots.savefig(pp,"out/r2pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon_hh_2.png")
        ind+= 1
    end
end

fctmain(32,256)
