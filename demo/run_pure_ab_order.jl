
include("../src/interface.jl")
using LinearAlgebra
using Plots
using Random

function fctmain(n_tau, prec)
    setprecision(prec)
    Random.seed!(8900161)
    u0 = rand(BigFloat, 4)
    B = 2rand(BigFloat, 4, 4)-ones(BigFloat,4, 4)
    alpha = 2rand(BigFloat, 4)-ones(BigFloat, 4)
    beta = 2rand(BigFloat, 4)-ones(BigFloat, 4)
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
    nbmaxtest=12
    ordmax=17
    debord=3
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    ind=1
    for order=debord:pasord:ordmax
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, order-debord+1)
        resnorm=0
        resnormprec=1
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
            resnorm = norm(diff,Inf)
            y[indc,ind] = min(norm(diff,Inf), 1.1)
            println("result=$y")
            println("log2(y)=$(log2.(y))")
            nb *= 2
            indc += 1
        end
        for i=debord:pasord:order
            labels[1,(i-debord)÷pasord+1] = " order=$i "
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
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res4_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact.pdf")
        if resnorm > resnormprec
            break
        end
        resnormprec = resnorm
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32, 512)
