include("../src/interface.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
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
    return ret
end
function fctmain(n_tau, prec)
    setprecision(prec)
    u0=BigFloat([90, -44, 83, 13]//100)
    t_max = big"1.0"
    epsilon=big"0.0017"
    println("epsilon=$epsilon")
    nbmaxtest=12
    ordmax=17
    debord=3
    pasord=2
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )

    res_gen = Array{ Array{BigFloat,1}, 2}(undef, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    prob = HiOscODEProblem(henon_heiles, u0, (big"0.0",t_max), missing, A, epsilon)
    tabsol = Array{Array{BigFloat,1},1}(undef,0)
    println("ordmax=$ordmax")
    indref = 1
    ind=1
    for order=debord:pasord:ordmax
        println("eps=$epsilon order=$order")
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, div(order-debord, pasord)+1)
        resnorm=0
        resnormprec=1
        sol =undef
        println("preparation ordre $order + 2")
        par_u0 = missing     
        while indc <= nbmaxtest
            res = solve(prob, nb_tau=n_tau, order=order, order_prep=ordprep, nb_t=nb,par_u0=par_u0, dense=false)
            par_u0=res.par_u0
            sol = res[end]          
            push!(tabsol, sol)
            res_gen[indc,ind] = sol 
            (a, b) = getmindif(tabsol)
            if a != 0
                indref = a
                solref = (tabsol[a]+tabsol[b])/2
                for i=1:ind
                    borne = (i <ind) ? size(res_gen,1) : indc
                    for j = 1:borne
                        nm = norm(res_gen[j,i] - solref, Inf)
                        y[j,i] = (nm < 1) ? nm : NaN
                   end
                end
            end
            x[indc] = 1.0/nb
            println("result=$y")
            println("log2(y)=$(log2.(y))")
            nb *= 2
            indc += 1
        end
        for i=debord:pasord:order
            labels[1,(i-debord)÷pasord+1] = " order = $i "
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
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res3p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.pdf")
        Plots.savefig(p,"out/res3p_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.png")
        Plots.savefig(pp,"out/res3pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.pdf")
        Plots.savefig(pp,"out/res3pp_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.png")
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32,512)
