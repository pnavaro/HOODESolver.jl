
include("../src/twoscales_pure_ab.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
function juliastd_solve(fct::Function, u0, t, eps, atol, rtol )
    prob = ODEProblem(henon_heiles,u0,(big"0.0",t),eps)
    sol = solve(prob, abstol=atol, reltol=rtol)
    return sol.u[end]
end

function twoscales_solve( par_u0::PrepareU0, order, t, nb)
    
    pargen = PrepareTwoScalePureAB(nb, t, order, par_u0)

    return twoscales_pure_ab(pargen, only_end=true)

end

function fctmain(n_tau)
    u0 =[big"0.1", big"0.11", big"0.15", big"0.10781"]
    B = [ big"-0.12984599677" big"-0.9277" big"0.32984110099677" big"0.142984599677"
    big"-0.4294599677" big"0.1273371193457" big"0.429841100997777222" big"0.99484599677"
    big"0.22984996779898" big"0.327667214311" big"0.1298410099677" big"-0.342984599677"
    big"0.7298459677881111007" big"-0.0278879898" big"0.7294110099677" big"-0.66294599677"
    ]
    t_max = big"1.0"
    epsilon=big"1.0"/big"12800.0"
    nbmaxtest=9
    ordmax=15
    debord=3
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)

    ind=1
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for order=debord:pasord:ordmax
        if order > 10
            setprecision(512)
        end
        u0=BigFloat.(u0)
        t_max = BigFloat(t_max)
        epsilon = BigFloat(epsilon)
        B = BigFloat.(B)
        fct = u -> B*u
        parphi = PreparePhi(epsilon, n_tau, A, fct, B)
        solref = getexactsol(parphi, u0, t_max)
        println("eps=$epsilon solRef=$solref order=$order")
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, order-debord+1)
        resnorm=0
        resnormprec=1
        println("preparation ordre $order + 2")
        @time par_u0 = PrepareU0(parphi, order+2, u0)       
        while indc <= nbmaxtest
            @time sol = twoscales_solve( par_u0, order, t_max, nb)
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
            labels[1,(i-debord)÷pasord+1] = " eps,order=$(convert(Float32,epsilon)),$i "
        end
        p=Plots.plot(
                        x,
                        view(y,:,1:ind),
                        xlabel="delta t",
                        xaxis=:log,
                        ylabel="error",
                        yaxis=:log,
                        legend=:bottomright,
                        label=labels,
                        marker=2
                    )
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/rapres_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_exact_v2g.pdf")
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
fctmain(32)
