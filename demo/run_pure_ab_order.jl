
include("../src/interface.jl")
using LinearAlgebra
using Plots
using Plots.PlotMeasures
using Random

function fctmain(n_tau, prec)
    setprecision(prec)
    Random.seed!(8900161)
    u0 = [-big"0.34",big"0.78",big"0.67",-big"0.56"]
    B = [big"0.12" -big"0.78" big"0.91" big"0.34"
    	-big"0.45" big"0.56" big"0.3" big"0.54"
    	-big"0.67" big"0.09" big"0.18" big"0.89"
    	-big"0.91" -big"0.56" big"0.11" -big"0.56"]


    alpha =  [big"0.12",-big"0.98",big"0.45",big"0.26"]
    beta =  [-big"0.4",big"0.48",big"0.23",-big"0.87"]
    u0=BigFloat.(u0)
    B = BigFloat.(B)
    println("u0=$u0")
    println("B=$B")
    println("alpha=$alpha")
    println("beta=$beta")
    fct = (u,p,t) -> B*u + t*p[1] +p[2]
    t_max = big"1.0"
    epsilon=big"0.015"
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    prob = HiOscODEProblem(fct,u0, (big"0.0",t_max), (alpha, beta), A, epsilon, B)
#    nbmaxtest=12
    nbmaxtest=7
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
#                        bottom_margin=50px
                   )
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res8essai_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.pdf")        
        Plots.savefig(p,"out/res8essai_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.png")
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32, 512)
