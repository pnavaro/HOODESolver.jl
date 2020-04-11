include("../src/interface.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
using Random

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
function fct(u, p, t)
    s1, c1 = sincos(u[1]/2)
    s2, c2 = sincos(u[2])
    s3, c3 = sincos(u[3])
    return [0, 0, u[6], c1*s2*s3/2, s1*c2*s3, s1*s2*c3]
end


function fctmain(n_tau, prec)

    setprecision(prec)
    Random.seed!(5612)
    t_max = big"1.0"
    epsilon=big"1e-7"
    println("epsilon=$epsilon")
    nbmaxtest=10
    ordmax=13
    debord=3
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    res_gen = Array{ Array{BigFloat,1}, 2}(undef, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    A = [0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 -1 0 0; 0 0 0 0 0 0]
    u0 = 2rand(BigFloat,6)-ones(BigFloat,6)
    t_0=big"0.0"
    t_max=big"1.0"
    prob = HiOscODEProblem(fct, u0, (t_0, t_max), missing, A, epsilon)
    nm = NaN
    ordmax += 1
    ordprep=undef
    nb = 20*2^(nbmaxtest)
    solref=undef
    while isnan(nm)
        ordmax -= 1
        @time sol = solve(prob, nb_t=nb, order=ordmax, getprecision=false, nb_tau=n_tau)
        solref = sol[end]
	    nm = norm(solref, Inf)
    end

    tabsol = Array{Array{BigFloat,1},1}(undef,1)

    tabsol[1] = solref

    println("ordmax=$ordmax")

    println("solref=$(solref[end])")

    solref_gen = solref
    indref = 1

    ind=1
    for order=debord:pasord:ordmax
        ordprep=order+2
        println("eps=$epsilon solRef=$(solref[end]) order=$order")
        nb = 10
        indc = 1
        labels=Array{String,2}(undef, 1, order-debord+1)
        resnorm=0
        resnormprec=1
        sol =undef
        println("preparation ordre $order + 2")
        while indc <= nbmaxtest
            @time solall = solve(prob, nb_t=nb, order=order, getprecision=false, nb_tau=n_tau)
            sol = solall[end]
            push!(tabsol, sol)
            res_gen[indc,ind] = sol
            (a, b), nm = getmindif(tabsol)
            if a != indref
                println("New solref !!!! a=$a, b=$b nm=$nm")
                indref = a
                solref = tabsol[a]
                for i=1:ind
                    borne = (i <ind) ? size(res_gen,1) : indc
                    for j = 1:borne
                        nm2 = min( norm(res_gen[j,i] - solref, Inf), 1.1)
                        y[j,i] = nm2 == 0 ? nm : nm2
                   end
                end
            else
                diff=solref-sol
                y[indc,ind] = min(norm(diff,Inf), 1.1)
            end
            x[indc] = 1.0/nb
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
                        legend=:topleft,
                        label=labels,
                        marker=2
                    )
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/r1_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_particle.pdf")
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
fctmain(32,512)
