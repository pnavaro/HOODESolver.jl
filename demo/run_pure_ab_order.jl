
include("../src/twoscales_pure_ab.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
using Random

function fctmain(n_tau, prec)
    setprecision(prec)
    Random.seed!(89161)
    u0 = rand(BigFloat, 4)
    B = 2rand(BigFloat, 4, 4)-ones(BigFloat,4, 4)
    t_max = big"1.0"
    epsilon=big"0.0005"
    nbmaxtest=11
    ordmax=17
    debord=3
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    ind=1
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    for order=debord:pasord:ordmax
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
        ordprep = min(order+2,10)
        ordprep = order+2
        println("preparation ordre $ordprep")
        @time par_u0 = PrepareU0(parphi, ordprep, u0)       
        while indc <= nbmaxtest
            @time pargen = PrepareTwoScalePureAB(nb, t_max, ordprep, par_u0)
            @time result, tabdf, tabdf2, nm, nm2 = twoscales_pure_ab(
    pargen,
    only_end=false,
    diff_fft=true
)
            pasaff = nb/100
            sol = result[:, end]
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            for i=1:50
                diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
                println("i=$i/$nb diff=$diff")
             end
             for i=51:pasaff:(nb-50)
                diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
                println("i=$i/$nb diff=$diff")
            end
            for i=(nb-50):nb
                diff = norm(result[:,i]-getexactsol(parphi,u0,t_max*(i-1)/nb))
                println("i=$i/$nb diff=$diff")
            end
            for i=1:50
                println("i=$i/$nb difffftInf=$(tabdf[i])")
            end
            pasaff=div(nb,100)
            for i=51:pasaff:(nb-50)
                println("i=$i/$nb difffftInf=$(tabdf[i])")
            end
            for i=(nb-50):nb
                println("i=$i/$nb difffftInf=$(tabdf[i])")
            end
            for i=1:50
                println("i=$i/$nb difffft2=$(tabdf2[i])")
            end
            for i=51:pasaff:(nb-50)
                println("i=$i/$nb difffft2=$(tabdf2[i])")
            end
            for i=(nb-50):nb
                println("i=$i/$nb difffft2=$(tabdf2[i])")
            end
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
        Plots.savefig(p,"out/res_$(prec_v)_$(eps_v)_$(order)_$(ordprep)_$(n_tau)_exact_v2h.pdf")
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
fctmain(32, 256)
