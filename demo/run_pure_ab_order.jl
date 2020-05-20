
include("../src/interface.jl")
using LinearAlgebra
using Plots
using Plots.PlotMeasures

function fctmain(n_tau, prec)
    setprecision(prec)
    u0 = BigFloat.([-34//100, 78//100, 67//100, -56//10])
    B = BigFloat.([12//100 -78//100 91//100 34//100
        -45//100 56//100 3//100 54//100
        -67//100 09//100 18//100 89//100
        -91//100 -56//100 11//100 -56//100])
    alpha =  BigFloat.([12//100, -98//100, 45//100, 26//100])
    beta =  BigFloat.([-4//100, 48//100, 23//100, -87//100])
    println("u0=$u0")
    println("B=$B")
    println("alpha=$alpha")
    println("beta=$beta")
    fct = (u,p,t) -> B*u + t*p[1] +p[2]
    t_max = big"1.0"
    epsilon=big"0.015"
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    prob = HiOscODEProblem(fct,u0, (big"0.0",t_max), (alpha, beta), A, epsilon, B)
#    nbmaxtest=6
    nbmaxtest=12
    ordmax=17
    debord=3
    pasord=2
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    ind=1
    for order=debord:pasord:ordmax
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, div(order-debord,pasord)+1)
        ordprep = min(order+2,10)
        ordprep = order+2
        println("preparation ordre $ordprep")
        par_u0=missing
        while indc <= nbmaxtest
            println("u0=$u0")
            println("B=$B")
            println("alpha=$alpha")
            println("beta=$beta")        
            @time res = solve(prob, nb_tau=n_tau, order=order, order_prep=ordprep, nb_t=nb,par_u0=par_u0, dense=false)
            par_u0=res.par_u0
            sol = res[end]
            solref=getexactsol(res.par_u0.parphi, u0, t_max)
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            diff=solref-sol
            x[indc] = 1.0/nb
            println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
            resnorm = Float64(norm(diff,Inf))
            y[indc,ind] = (resnorm < 1) ? resnorm : NaN
            println("result=$y")
            println("log2(y)=$(log2.(y))")
            nb *= 2
            indc += 1
        end
        for i=debord:pasord:order
            labels[1,div((i-debord),pasord)+1] = " order=$i "
        end
        yv = y[:,1:ind]
        p=Plots.plot(
    x,
    yv,
    xlabel="Δt",
    xaxis=:log,
    ylabel="error",
    yaxis=:log,
    legend=:bottomright,
    label=labels,
    marker=2,
)
        pp=Plots.plot(
    x,
    yv,
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
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res9p_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.pdf")        
        Plots.savefig(p,"out/res9p_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.png")
        Plots.savefig(pp,"out/res9pp_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.pdf")        
        Plots.savefig(pp,"out/res9pp_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.png")
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32, 512)
