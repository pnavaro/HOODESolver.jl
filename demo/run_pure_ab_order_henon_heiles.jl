include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Plots
using Random


# function getindmin( tab::Array{Array{BigFloat,1},1} )

    # summin=Inf
    # for i=1:size(tab,1)
        # sum = zero(BigFloat)
        # for j=i:size(tab,1)
            # sum += norm(tab[i]-tab[j],Inf)
        # end
        # if sum < summin
            # sum = summin
            # ret = i
        # end
    # end
# end



# function getsolref( solref_gen, res_gen, y)
    # mindiff = 1.0
    # for i=1:size(y,1), j=1:size(y,2)
        # if y[i,j] != NaN && y[i,j] != undef
            # diff = norm(solref_gen-res_gen[i,j], Inf)
            # if diff < mindiff && diff != 0.0
                # mindiff = diff
# 
		# end
# end
# end
# end

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


function fctmain(n_tau, prec)

    setprecision(prec)
    Random.seed!(5612)
    u0=rand(BigFloat,4)

    t_max = big"1.0"
    epsilon=big"0.001"
    println("epsilon=$epsilon")
    nbmaxtest=8
    ordmax=10
    debord=3
    pasord=1
    y = ones(Float64, nbmaxtest, div(ordmax-debord,pasord)+1 )

    res_gen = Array{ Array{BigFloat,1}, 2}(undef, nbmaxtest, div(ordmax-debord,pasord)+1 )
    x=zeros(Float64,nbmaxtest)
    A = [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    parphi = PreparePhi(epsilon, n_tau, A, henon_heiles)

    nm = NaN

    ordmax += 1
    ordprep=undef
    nb = 100*2^(nbmaxtest)
    solref=undef

    while isnan(nm)

        ordmax -= 1
        ordprep = ordmax+2

    	@time par_u0 = PrepareU0(parphi, ordprep, u0)
        @time pargen = PrepareTwoScalePureAB(nb, t_max, ordmax, par_u0)
        @time solref= twoscales_pure_ab(
pargen,
only_end=true,
diff_fft=false
)
	    nm = norm(solref, Inf)
    end

    tabsol = Array{Array{BigFloat,1},1}(undef,1)

    tabsol[1] = solref

    println("ordmax=$ordmax")

    println("solref=$solref")
    solref_gen = solref
    indref = 1

    ind=1
    for order=debord:pasord:ordmax
        ordprep=order+2
        println("eps=$epsilon solRef=$solref order=$order")
        nb = 100
        indc = 1
        labels=Array{String,2}(undef, 1, order-debord+1)
        resnorm=0
        resnormprec=1
        sol =undef
        println("preparation ordre $order + 2")
        @time par_u0 = PrepareU0(parphi, order+2, u0)       
        while indc <= nbmaxtest
            @time pargen = PrepareTwoScalePureAB(nb, t_max, ordprep, par_u0)
            @time sol = twoscales_pure_ab(
    pargen,
    only_end=true,
    diff_fft=false
)
            push!(tabsol, sol)
            res_gen[indc,ind] = sol 
            (a, b), nm = getmindif(tabsol)
            if a != indref
                println("New solref !!!! a=$a, b=$b nm=$nm")
                indref = a
                solref = tabsol[a]
                for i=1:indc
                    borne = (i <indc) ? size(res_gen,1) : ind
                    for j = 1:borne
                        nm2 = min( norm(res_gen[i,j] - solref, Inf), 1.1)
                        y[i,j] = nm2 == 0 ? nm : nm2
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
                        legend=:bottomright,
                        label=labels,
                        marker=2
                    )
        prec_v = precision(BigFloat)
        eps_v = convert(Float32,epsilon)
        Plots.savefig(p,"out/res2_$(prec_v)_$(eps_v)_$(order)_$(n_tau)_henon_heiles.pdf")
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
fctmain(64,512)
