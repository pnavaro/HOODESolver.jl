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

function getmindif(tab::Vector{Vector{BigFloat}})
    nmmin = Inf
    ret = (1, 2)
    for i=1:size(tab,1)
        for j=(i+1):size(tab,1)
            nm = norm(tab[i]-tab[j], Inf)
            if nm < nmmin
                nmmin = nm
                ret = i, j
            end
        end
    end
    return ret, min(nmmin,1.1)
end

function fctMain(n_tau)

 #   u0 =[big"0.12345678", big"0.1209182736", big"0.1290582671", big"0.1239681094" ]
    seed=12996
    Random.seed!(seed)
    u0=rand(BigFloat,4)
    println("seed = $seed")
 #   tab_eps = zeros(BigFloat,7)
 #   tab_eps= [big"0.5", big"0.2", big"0.1",big"0.05", big"0.02", big"0.01",big"0.005", big"0.002", big"0.001",big"0.0005", big"0.0002", big"0.0001"]
    tab_eps= [big"0.1",big"0.01", big"0.001",big"0.0001", big"0.00001", big"0.000001", big"0.0000001", big"0.00000001", big"0.000000001"]
    # epsilon=big"0.1"
    # for i=1:7
    #     tab_eps[i] = epsilon
    #     epsilon /= 100
    # end
    nbmaxtest=11
    t_max = big"1.0"
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    tabtabsol = Vector{Vector{Vector{BigFloat}}}()
    for order=2:14
        ordprep=order+2
        y = ones(Float64, nbmaxtest, size(tab_eps,1) )
        y_big = ones(BigFloat, nbmaxtest, size(tab_eps,1) )
        x=zeros(Float64,nbmaxtest)
        ind=1    
        for epsilon in tab_eps
    #        fct = u -> [ 1.11u[2]^2+0.8247u[4]^3+0.12647u[3]^4, 
    #       0.4356789141u[1]-0.87u[3]*u[4], 
    #       1.9898985/(1+1.1237u[4]^3), 0.97u[1]*u[3]+0.8111u[2]^2 ]
            fct = (u,p,t) ->  [ 
    u[3]^3 + sin(t^2), 
    u[4] + u[1]^2 + 1/(1+exp(t))-0.35, 
    u[1]*u[2]*u[4]*(1+t), 
    u[2] - u[2]^2 + 1/(1+u[4]^2) + u[1]^2 + cos(exp(t))
]
            parphi = PreparePhi(epsilon, n_tau, A, fct)
            println("prepareU0 eps=$epsilon n_tau=$n_tau")
            nb=10
            @time par_u0 = PrepareU0(parphi, ordprep, u0)
            @time pargen = PrepareTwoScalesPureAB(nb*2^nbmaxtest, t_max, order, par_u0)
            @time solref = twoscales_pure_ab(pargen, only_end=true)
            if size(tabtabsol,1) < ind
                push!(tabtabsol,Vector{Vector{BigFloat}}())
            end
            tabsol = tabtabsol[ind]
            push!(tabsol, solref)
            res_gen = Array{ Array{BigFloat,1}, 1}(undef, nbmaxtest)
            indref = 1   
            eps_v = convert(Float32,epsilon)
            println("epsilon = $eps_v solref=$solref")
            indc =1
            labels=Array{String,2}(undef, 1, ind)  
            while indc <= nbmaxtest
                @time pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
                @time sol= twoscales_pure_ab(pargen, only_end=true)
                push!(tabsol, sol)
                res_gen[indc] = sol
                diff=solref-sol
                println("tabsol=$tabsol")
                (a, b), nm = getmindif(tabsol)
                if a != indref
                    println("New solref !!!! a=$a, b=$b nm=$nm")
                    indref = a
                    solref = tabsol[a]
                    for i=1:indc
                        nm2 = min( norm(res_gen[i] - solref, Inf), 1.1)
                        y[i, ind] = nm2 == 0 ? nm : nm2
                        y_big[i, ind] = nm2 == 0 ? nm : nm2
                    end
                else
                    diff=solref-sol
                    y[indc,ind] = min(norm(diff,Inf), 1.1)
                    y_big[indc,ind] = min(norm(diff,Inf), 1.1)
                end
                println("solref=$solref")
                println("nb=$nb sol=$sol")
                x[indc] = 1.0/nb
                println("epsilon=$epsilon result=$y")
                println("epsilon=$epsilon reslog2=$(log2.(y))")
                println("epsilon=$epsilon result=$y_big")
                println("epsilon=$epsilon reslog2=$(log2.(y_big))")
                nb *= 2
                indc += 1
            end
            if ind == size(tab_eps,1)
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
                Plots.savefig(p,"out/r20_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_epsilon_fct.pdf")
            end
            ind+= 1
        end
    end
end


# testODESolverEps()

# for i=3:9
#     fctMain(2^i)
# end
# setprecision(512)
fctMain(32)
