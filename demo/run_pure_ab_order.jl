
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
    u0 =[big"1"/8, big"1"/8+big"1"/128, big"1"/8+big"2"/128, big"1"/8+big"3"/128]
    epsilon=big"1"/256
    nbmaxtest=9
    ordmax=5
    debord=3
    y = ones(Float64, nbmaxtest, div(ordmax-debord,2)+1 )
 #   y = ones(Float64, nbmaxtest, size(tabEps,1) )
    x=zeros(Float64,nbmaxtest)

    ind=1
 
    parphi = PreparePhi(epsilon, n_tau, [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], henon_heiles)
    par_u0 = PrepareU0(parphi, ordmax+2, u0)    
	nball = 100*2^(nbmaxtest+1)
    @time solrefall = twoscales_solve( par_u0, ordmax, big"1.0", nball)
	solref = solrefall[:,end]
    for order=debord:ordmax
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
        for i=debord:2:order
            labels[1,(i-debord)÷2+1] = " eps,order=$(convert(Float32,epsilon)),$i "
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
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_all+v5.pdf")
        ind+= 1
    end
end

# testODESolver()

# for i=3:9
#     fctMain(2^i)
# end
fctmain(32)
