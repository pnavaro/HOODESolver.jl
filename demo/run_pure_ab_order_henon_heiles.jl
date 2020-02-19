
include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
using Random
function twoscales_solve( par_u0::PrepareU0, order, t, nb)
    
    pargen = PrepareTwoScalePureAB(nb, t, order, par_u0)

    return twoscales_pure_ab(pargen, only_end=true)

end

function fctmain(n_tau)

    Random.seed!(5612)
    u0=rand(BigFloat,4)

    t_max = big"1.0"
    epsilon=big"0.00005"
    nbmaxtest=7
    ordmax=12
    debord=2
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    parphi = PreparePhi(epsilon, n_tau, A, henon_heiles)
    @time par_u0 = PrepareU0(parphi, ordmax+2, u0)
    
    @time solref = twoscales_solve( par_u0, ordmax, t_max, 100*2^(nbmaxtest+1))

    println("solref=$solref")

    
    

    ind=1
    for order=debord:pasord:ordmax
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
        Plots.savefig(p,"out/res_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.pdf")
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
