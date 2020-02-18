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
    
    parGen = PrepareTwoScalePureAB(nb, t, order, par_u0)

    return twoscales_pure_ab(parGen, only_end=true)

end


function fctMain(n_tau)

    u0 =[big"0.12345678", big"0.1209182736", big"0.1290582671", big"0.1239681094" ]

    Random.seed!(123456789)

    B = 2rand(BigFloat, 4, 4) - ones(BigFloat, 4, 4)
    println("B=$B")

    tab_eps = zeros(BigFloat,8)
    epsilon=big"0.1"
    for i=1:8
        tab_eps[i] = epsilon
        eps /= 10
    end

#    tab_eps=[big"0.005"]
    nbMaxTest=8
    order=6
    t_max = big"1.0"
 #    y = ones(Float64, nbMaxTest, ordmax-debord+1 )
    y = ones(Float64, nbMaxTest, size(tab_eps,1) )
    x=zeros(Float64,nbMaxTest)

    ind=1
    for epsilon in tab_eps
  #      @time solref = juliaStdSolve(henonHeiles, u0, big"1.0", eps, 1e-40, 1e-40 )

        fct = u -> B*u
        parphi = PreparePhi(n_tau, epsilon, [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], fct, B)
        println("prepareU0 eps=$epsilon n_tau=$n_tau")
        @time par_u0 = PrepareU0(parphi, 2, u0)

        solref = getexactsol(parphi, u0, t_max)

        eps_v = convert(Float32,epsilon)

        println("epsilon = $eps_v solref=$solref")
        nb = 100
        indc =1
        labels=Array{String,2}(undef, 1, ind)  
        while indc <= nbMaxTest
            @time sol = twoscales_solve( par_u0, order, big"1.0", nb) 
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
    legend=:bottomright,
    label=labels,
    marker=2
)
        
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,eps)
        Plots.savefig(p,"out/resAB_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon_v1.pdf")
        ind+= 1
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctMain(32)
