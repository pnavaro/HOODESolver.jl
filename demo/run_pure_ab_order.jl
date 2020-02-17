
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

    return twoscales_pure_ab(pargen)

end

function fctmain(n_tau)
    u0 =[big"0.1", big"0.11", big"0.15", big"0.10781"]
    B = [ big"-0.12984599677" big"-0.9277" big"0.32984110099677" big"0.142984599677"
    big"-0.4294599677" big"0.127337" big"0.4298411009977" big"0.99484599677"
    big"0.2298499677" big"0.327667" big"0.1298410099677" big"-0.342984599677"
    big"0.7298459677" big"-0.027887" big"0.7294110099677" big"-0.66294599677"
    ]

    epsilon=big"0.001"
    nbmaxtest=11
    ordmax=16
    debord=10
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )
 #   y = ones(Float64, nbmaxtest, size(tabEps,1) )
    x=zeros(Float64,nbmaxtest)

    ind=1
    fct = u -> B*u
    parphi = PreparePhi(epsilon, n_tau, [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], fct, B)
    println("Preparation ordre 2")
    @time par_u0 = PrepareU0(parphi, 2, u0)
    println("fin preparation")
	nball = 100*2^nbmaxtest
    solrefall = zeros(4,nball+1)
    for i=1:nball
        solrefall[:,i+1] = getexactsol(parphi, u0, i*big"1.0"/nball)
    end
    solrefall[:,1] = u0
    solref = solrefall[:,end]
    for order=debord:pasord:ordmax
        println("eps=$epsilon solRef=$solref order=$order")
        nb = 100
        indc =1
        labels=Array{String,2}(undef, 1, order-debord+1)
        par_u0 = PrepareU0(parphi, order+2, u0)       
        while indc <= nbmaxtest
            @time resall = twoscales_solve( par_u0, order, big"1.0", nb)
            coef = div(nball,nb)
            for i=1:(nb+1)
                println("tr__$nb i=$i error=$(norm(resall[:,i]-solrefall[:,coef*(i-1)+1]))")
            end
            sol = resall[:,end]
            println("solref=$solref")
            println("nb=$nb sol=$sol")
            diff=solref-sol
            x[indc] = 1.0/nb
            println("nb=$nb dt=$(1.0/nb) normInf=$(norm(diff,Inf)) norm2=$(norm(diff))")
            y[indc,ind] = norm(diff,Inf)
            println("result=$y")
            nb *= 2
            indc += 1
        end
        for i=debord:pasord:order
            labels[1,(i-debord)÷pasord+1] = " eps,order=$(convert(Float32,epsilon)),$i "
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
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_exact.pdf")
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32)
