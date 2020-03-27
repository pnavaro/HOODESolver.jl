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

 #   u0 =[big"0.12345678", big"0.1209182736", big"0.1290582671", big"0.1239681094" ]
    seed=1293456
    Random.seed!(seed)
    tab_eps = zeros(BigFloat,6)
    epsilon=big"0.1"
    for i=1:6
        tab_eps[i] = epsilon
        epsilon /= 10
    end
    nbMaxTest=9
    order=7
    ordprep=order+2
    t_max = big"3.0"
    y = ones(Float64, nbMaxTest, size(tab_eps,1) )
    x=zeros(Float64,nbMaxTest)
    ind=1
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for epsilon in tab_eps
        u0=rand(BigFloat,4)
        B = 2rand(BigFloat, 4, 4) - ones(BigFloat, 4, 4)
        paramfct = (2rand(BigFloat, 4) - ones(BigFloat, 4),
                        2rand(BigFloat, 4) - ones(BigFloat, 4))
        fct = (u, p, t) -> B*u + t*p[1]+p[2]
        parphi = PreparePhi(epsilon, n_tau, A , fct, B, paramfct=paramfct )
#       fct = (u, p, t) -> B*u
#       parphi = PreparePhi(epsilon, n_tau, A , fct, B)
        println("prepareU0 eps=$epsilon n_tau=$n_tau")
        @time par_u0 = PrepareU0(parphi, ordprep, u0)
        solref = getexactsol(parphi, u0, t_max)
        eps_v = convert(Float32,epsilon)
        println("epsilon = $eps_v solref=$solref")
        nb = 30
        indc =1
        labels=Array{String,2}(undef, 1, ind)  
        while indc <= nbMaxTest
            @time pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
            @time result, tabdf, tabdf2, nm, nm2 = twoscales_pure_ab(
    pargen,
    only_end=false,
    diff_fft=true
)
            sol = result[:, end]
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            # pasaff=div(nb,100)
            # for i=1:50
            #     diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
            #     println("i=$i/$nb diff=$diff")
            #  end
            #  for i=51:pasaff:(nb-50)
            #     diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
            #     println("i=$i/$nb diff=$diff")
            # end
            # for i=(nb-50):nb
            #     diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
            #     println("i=$i/$nb diff=$diff")
            # end
            # for i=1:50
            #     println("i=$i/$nb difffftInf=$(tabdf[i])")
            # end
            # for i=51:pasaff:(nb-50)
            #     println("i=$i/$nb difffftInf=$(tabdf[i])")
            # end
            # for i=(nb-50):nb
            #     println("i=$i/$nb difffftInf=$(tabdf[i])")
            # end
            # for i=1:50
            #     println("i=$i/$nb difffft2=$(tabdf2[i])")
            # end
            # for i=51:pasaff:(nb-50)
            #     println("i=$i/$nb difffft2=$(tabdf2[i])")
            # end
            # for i=(nb-50):nb
            #     println("i=$i/$nb difffft2=$(tabdf2[i])")
            # end
            diff=solref-sol
            x[indc] = t_max/nb
            println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
            y[indc,ind] = min(norm(diff,Inf),big"1.0")
            println("epsilon=$epsilon result=$y")
            println("epsilon=$epsilon reslog2=$(log2.(y))")
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
        Plots.savefig(p,"out/r4_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_epsilon.pdf")
        ind+= 1
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctMain(256)
