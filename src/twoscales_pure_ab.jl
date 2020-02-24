include("preparephi.jl")
include("coefexp_ab.jl")

"""
    PrepareTwoScalePureAB(n_max, t_max, order, par_u0::PrepareU0)


"""


struct PrepareTwoScalePureAB
    n_max
    t_max
    dt
    order
    parphi::PreparePhi
    par_u0::PrepareU0
    p_coef::CoefExpAB
    exptau
    exptau_inv   
    function PrepareTwoScalePureAB(n_max, t_max, order, par_u0::PrepareU0)
        parphi = par_u0.parphi
        T = typeof(parphi.epsilon)
        dt = T(t_max)/n_max
        p_coef = CoefExpAB(order, parphi.epsilon, parphi.n_tau, dt)
        exptau = collect(transpose(exp.(-im*dt / parphi.epsilon * parphi.tau_list)))
        exptau_inv = collect(transpose(exp.(im*dt / parphi.epsilon * parphi.tau_list)))
        return new(
    n_max, 
    t_max, 
    dt, 
    order, 
    parphi,
    par_u0,
    p_coef, 
    exptau,
    exptau_inv
)
    end
end

function _calculfft(par::PrepareTwoScalePureAB, resfft)
    f = filtredfct(par.parphi, real(ifftgen(par.parphi.par_fft, resfft)))
    return fftgen(par.parphi.par_fft, f)
end
solJul = zeros(BigFloat,4)
solJul_neg = zeros(BigFloat,4)

function _calcul_ab(par::PrepareTwoScalePureAB, ord, fftfct, fft_u, dec, sens)
 # println("_calcul_ab par.order=$(par.order) ord=$ord dec=$dec sens=$sens dt=$(par.dt) xxxx")
    resfft = fft_u[dec-sens] .* ((sens==1) ? par.exptau : par.exptau_inv)
    tab_coef = (sens == 1) ? par.p_coef.tab_coef : par.p_coef.tab_coef_neg
    for k=1:ord
        indice = dec-sens*k
#        println("k=$k indice=$indice")
# println("size(tab_coef,1)=$(size(tab_coef,1)) size(fftfct[indice])=$(size(fftfct[indice]))")
        resfft .+= transpose(tab_coef[:, k, ord]).*fftfct[indice]
    end
    fft_u[dec] = resfft
    fftfct[dec] = _calculfft(par, resfft)

    # if isexactsol(par.parphi)
    #     u = _getresult(fft_u[dec], (dec-par.order)*par.dt, par.parphi)
    #     u_exact = getexactsol(par.parphi, par.par_u0.u0, (dec-par.order)*par.dt)
    #     println("ordre = $ord/$(par.order) indice=$(dec-par.order) normdiff=$(norm(u-u_exact,Inf))")
    # end

    # if dec == par.order+sens
    #     u = _getresult(fft_u[par.order+sens], sens*par.dt, par.parphi)
        # println("u=$(convert(Array{Float64,1},u))")
        # if sens == 1
        #     global solJul
        #     println("solJul=$(convert(Array{Float64,1},solJul))")
        #     println("sens=$sens par.order=$(par.order) ord=$ord diff=$(norm(solJul-u))")
        # else
        #     global solJul_neg
        #     println("solJul_neg=$(convert(Array{Float64,1},solJul_neg))")
        #     println("sens=$sens par.order=$(par.order) ord=$ord diff=$(norm(solJul_neg-u))")
        # end
    # end
end

function _init_ab(par::PrepareTwoScalePureAB, fftfct, fft_u)
    println("_init_ab order=$(par.order)")
    if par.order != 1
#        oo=0
        for new_ord=2:par.order #  in vcat(collect(2:par.order),[par.order])
 #           if oo != new_ord
                _calcul_ab(par, new_ord-1, fftfct, fft_u, par.order+new_ord-1, 1)
 #           end
            for k = 1:new_ord-1
                _calcul_ab(par, new_ord, fftfct, fft_u, par.order-k, -1)
            end
            for k = 1:new_ord-1
                _calcul_ab(par, new_ord, fftfct, fft_u, par.order+k, 1)
            end
 #           oo=new_ord
        end
    end
end
# function tr_ab(vect, ut0, ord, par::PrepareTwoScaleGen)
 
#     resfft = par.exptau.*fftGen(par.parphi.par_fft, ut0)

#     f = filtredFctGen(ut0, par.parphi)
#     f = fftGen(par.parphi.par_fft, f)

#     resfft .+= par.precompute.tab_coef[1,ord,:].*f

#     for i=1:ord-1
#         resfft .+= par.precompute.tab_coef[i+1,ord,:].*vect[i]
#     end

#     res = real(ifftGen(par.parphi.par_fft, resfft))

#     return res, f
# end

function _tr_ab(par::PrepareTwoScalePureAB, fftfct, u_chap)
    resfft = par.exptau.* u_chap
    bound = par.order-1
    for k =0:bound
#        println("k=$k coef[1,2]=$(par.p_coef.tab_coef[k+1,par.order,1:2])")
        resfft .+= transpose(par.p_coef.tab_coef[:, k+1, par.order]).*fftfct[end-k]
    end
    f = _calculfft(par, resfft)
    return f, resfft
end
function _getresult(u_chap, t, par::PreparePhi)
    # matlab : u1=real(sum(fft(ut)/Ntau.*exp(1i*Ltau*T/ep)));
    u1 = real(u_chap * exp.(1im * par.tau_list * t / par.epsilon)) / par.n_tau
#    res = exp(getA(par.parPhi) * (n - 1) * par.dt / par.parPhi.eps) * u1
    res = exp( t / par.epsilon *par.sparse_A) * u1
    return res
end


function twoscales_pure_ab(par::PrepareTwoScalePureAB; only_end=false, borne_fft=1.0)

    fftfct = Vector{Array{Complex{BigFloat}, 2}}(undef, 2par.order-1)

    fft_u = Vector{Array{Complex{BigFloat}, 2}}(undef, 2par.order-1)

    
    res_u = par.par_u0.ut0
    u0 = par.par_u0.u0
   
    fftfct[par.order] = fftgen(par.parphi.par_fft, filtredfct(par.parphi, res_u))

    println("twoscales_pure_ab epsilon/dt=$(par.parphi.epsilon/par.dt)")
    println("twoscales_pure_ab dt/epsilon=$(par.dt/par.parphi.epsilon)")
    println("twoscales_pure_ab order=$(par.order)")

  #  dump(par)

  #  println("res_u=$res_u")
    fft_u[par.order] = fftgen(par.parphi.par_fft, res_u)

    # println("res[0]=$(convert(Array{Float64,1},_getresult(fft_u[par.order],0,par.parphi)))")
    # println("res[0]bis=$(convert(Array{Float64,1}, _getresult(real(ifftgen(par.parphi.par_fft, par.exptau.*fftgen(par.parphi.par_fft,fft_u[par.order]))),0,par.parphi)))")

    result = zeros(typeof(par.parphi.epsilon), par.parphi.size_vect, par.n_max+1)
    result[:,1] = u0

   # println("diff result ori = $(norm(u0-_getresult(fft_u[par.order], 0, par.parphi)))")

    _init_ab(par, fftfct, fft_u)

    memfft = fftfct[par.order:(2par.order-1)]

    if !only_end
        for i=2:par.order
            result[:,i] = _getresult(fft_u[par.order+i-1], (i-1)*par.dt, par.parphi)
        end
    end

     # ring permutation where the beginning becomes the end and the rest is shifted by one
    permut = collect(Iterators.flatten((2:par.order,1:1)))

    ut0_fft = fft_u[end]
    println("")
    norm_delta_fft = 0
    nbnan = 0
    borne_nm=0
    c_mult=1.1
    for i=par.order:par.n_max
        if i%10000 == 1
            println(" $(i-1)/$(par.n_max)")
        end
        resfft, ut0_fft = _tr_ab(par, memfft, ut0_fft)
        if !only_end
            result[:,i+1] = _getresult(ut0_fft, i*par.dt, par.parphi)
        end
        nm = norm(resfft-memfft[end],Inf)
        norm_delta_fft = max(nm,norm_delta_fft)
        if ( nm > borne_nm || isnan(nm) ) && nbnan < 10
            borne_nm = nm*c_mult
            c_mult *= 2
            println("i=$i nm=$nm")
            nbnan += isnan(nm) ? 1 : 0
        end
        memfft = memfft[permut]
        memfft[end] = resfft
        if i%100 == 0
            print("x")
        end
    end

    println("norm diff fft = $norm_delta_fft")

    return only_end ? _getresult(ut0_fft, par.t_max, par.parphi) : result
end
