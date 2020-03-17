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
    ret = (0, 0)
    for i=1:size(tab,1)
        for j=(i+1):size(tab,1)
            nm = norm(tab[i]-tab[j], Inf)
            if nm < nmmin
                nmmin = nm
                ret = i, j
            end
        end
    end
    return ret, nmmin
end

function fctMain(n_tau)

 #   u0 =[big"0.12345678", big"0.1209182736", big"0.1290582671", big"0.1239681094" ]
    seed=126
    Random.seed!(seed)
    u0=rand(BigFloat,4)
    println("seed = $seed")
    tab_eps = zeros(BigFloat,7)
    epsilon=big"0.1"
    for i=1:7
        tab_eps[i] = epsilon
        epsilon /= 100
    end
    nbmaxtest=8
    order=4
    t_max = big"1.0"
    A=[0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for order=5:17
        ordprep=order+2
        y = ones(Float64, nbmaxtest, size(tab_eps,1) )
        y_big = ones(BigFloat, nbmaxtest, size(tab_eps,1) )
        x=zeros(Float64,nbmaxtest)
        ind=1    
        for epsilon in tab_eps
    #        fct = u -> [ 1.11u[2]^2+0.8247u[4]^3+0.12647u[3]^4, 
    #       0.4356789141u[1]-0.87u[3]*u[4], 
    #       1.9898985/(1+1.1237u[4]^3), 0.97u[1]*u[3]+0.8111u[2]^2 ]
            fct = u -> [u[3]^3, u[4]+u[2]+u[3]^2.61, 
                        u[1]*u[2], u[2]-u[2]^2+1/(1+u[4]^2)+u[1]^1.77 ]
            parphi = PreparePhi(epsilon, n_tau, A, fct)
            println("prepareU0 eps=$epsilon n_tau=$n_tau")
            nb=10
            @time par_u0 = PrepareU0(parphi, ordprep, u0)
            @time pargen = PrepareTwoScalePureAB(nb*2^nbmaxtest, t_max, order, par_u0)
            @time solref = twoscales_pure_ab(pargen, only_end=true)
            tabsol = Array{Array{BigFloat,1},1}(undef,1)
            res_gen = Array{ Array{BigFloat,1}, 1}(undef, nbmaxtest)

            tabsol[1] = solref
            indref = 1   
            eps_v = convert(Float32,epsilon)
            println("epsilon = $eps_v solref=$solref")
            indc =1
            labels=Array{String,2}(undef, 1, ind)  
            while indc <= nbmaxtest
                @time pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
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
            legend=:topleft,
            label=labels,
            marker=2
        )
                
                prec_v = precision(BigFloat)
                Plots.savefig(p,"out/r5_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_epsilon_fct.pdf")
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
