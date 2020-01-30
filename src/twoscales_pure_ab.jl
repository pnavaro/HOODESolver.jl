include("fftBig.jl")
include("polylagrange.jl")
include("../src/henonHeiles.jl")
include("twoScale3.jl")

struct PrepareTwoScalePureAB
    n_max
    t_max
    dt
    order
    parphi::PreparePhi
    par_u0::PrepareU0
    p_lagr::PolyLagrange
    expTau
    expTauInv
    
    function PrepareTwoScalePureAB(n_max, t_max, order, par_u0::PrepareU0)
        parphi = par_u0.parphi
        dt = t_max/n_max
        N, coef = typeof(parphi.eps) == BigFloat ? (BigInt, 10) : (Int64, 1)
        eps_rat = rationalize(N, parphi.eps, coef*eps(parphi.eps))
        coef_rat = rationalize.(N, parphi.coefTauList)
        dt_rat = rationalize(N, dt, coef*eps(dt))

println("eps_rat=$eps_rat coef_rat=$coef_rat dt_rat=$dt_rat")
        p_lagr = PolyLagrange(order, eps_rat, coef_rat, dt_rat)
        expTau = exp.(-im*dt / parphi.eps * parphi.coefTauList)
        expTauInv = exp.(im*dt / parphi.eps * parphi.coefTauList)

        return new(
    n_max, 
    t_max, 
    dt, 
    order, 
    parphi,
    par_u0,
    p_lagr, 
    expTau,
    expTauInv
)
    end
end

function _calculfft(par::PrepareTwoScalePureAB, resfft)
    f = filtredFctGen(real(ifftGen(par.parphi.parFft, resfft)), par.parphi)
    return fftGen(par.parphi.parFft, f)
end
solJul = zeros(BigFloat,4)
solJul_neg = zeros(BigFloat,4)

function _calculAB(par::PrepareTwoScalePureAB, ord, fftfct, fft_u, dec, sens)
 println("_calculAB par.order=$(par.order) ord=$ord dec=$dec sens=$sens dt=$(par.dt) xxxx")
    resfft = ((sens==1) ? par.expTau : par.expTauInv).*fft_u[dec-sens]
    respe = (sens == 1) ? par.p_lagr.respe : par.p_lagr.respe_neg
    for k=1:ord
        indice = dec-sens*k
        println("k=$k indice=$indice")
        resfft .+= respe[k,ord,:].*fftfct[indice]
    end
    fft_u[dec] = resfft
    fftfct[dec] = _calculfft(par, resfft)

    if dec == par.order+sens
        u = _getresult(fft_u[par.order+sens], sens*par.dt, par.parphi)
        println("u=$(convert(Array{Float64,1},u))")
        if sens == 1
            global solJul
            println("solJul=$(convert(Array{Float64,1},solJul))")
            println("sens=$sens par.order=$(par.order) ord=$ord diff=$(norm(solJul-u))")
        else
            global solJul_neg
            println("solJul_neg=$(convert(Array{Float64,1},solJul_neg))")
            println("sens=$sens par.order=$(par.order) ord=$ord diff=$(norm(solJul_neg-u))")
        end
    end
end

function _initAB(par::PrepareTwoScalePureAB, fftfct, fft_u)
    println("_initAB order=$(par.order)")
    if par.order != 1
        for new_ord = 2:par.order
            _calculAB(par, new_ord-1, fftfct, fft_u, par.order+new_ord-1, 1)
            for k = 1:new_ord-1
                _calculAB(par, new_ord, fftfct, fft_u, par.order-k, -1)
            end
            for k = 1:new_ord-1
                _calculAB(par, new_ord, fftfct, fft_u, par.order+k, 1)
            end
        end
    end
end
# function trAB(vect, ut0, ord, par::PrepareTwoScaleGen)
 
#     resfft = par.expTau.*fftGen(par.parphi.parFft, ut0)

#     f = filtredFctGen(ut0, par.parphi)
#     f = fftGen(par.parphi.parFft, f)

#     resfft .+= par.precompute.respe[1,ord,:].*f

#     for i=1:ord-1
#         resfft .+= par.precompute.respe[i+1,ord,:].*vect[i]
#     end

#     res = real(ifftGen(par.parphi.parFft, resfft))

#     return res, f
# end

function _trAB(par::PrepareTwoScalePureAB, fftfct, u_chap)
    resfft = par.expTau.* u_chap
    bound = par.order-1
    for k =0:bound
#        println("k=$k coef[1,2]=$(par.p_lagr.respe[k+1,par.order,1:2])")
        resfft .+= par.p_lagr.respe[k+1,par.order,:].*fftfct[end-k]
    end
    f = _calculfft(par, resfft)
    return f, resfft
end
function _getresult(u_chap, t, par::PreparePhi)
    # matlab : u1=real(sum(fft(ut)/Ntau.*exp(1i*Ltau*T/ep)));
    u1 = real(sum(u_chap .* exp.(1im * par.coefTauList * t / par.eps);
                    dims = (1,)
                )) / par.nTau

    u1 = reshape(u1,(4,))

#    res = exp(getA(par.parPhi) * (n - 1) * par.dt / par.parPhi.eps) * u1
    res = expA(par, t / par.eps) * u1
    return res
end


function twoScalePureAB(par::PrepareTwoScalePureAB)

    fftfct = Vector{Union{Nothing, Array{Complex{BigFloat}, 2}}}(nothing, 2par.order-1)

    fft_u = Vector{Union{Nothing, Array{Complex{BigFloat}, 2}}}(nothing, 2par.order-1)

    
    res_u = par.par_u0.ut0
    u0 = par.par_u0.u0
   
    fftfct[par.order] = fftGen(par.parphi.parFft, filtredFctGen(res_u, par.parphi))

  #  dump(par)

  #  println("res_u=$res_u")
    fft_u[par.order] = fftGen(par.parphi.parFft, res_u)

    println("res[0]=$(convert(Array{Float64,1},_getresult(fft_u[par.order],0,par.parphi)))")
    println("res[0]bis=$(convert(Array{Float64,1}, _getresult(real(ifftGen(par.parphi.parFft, par.expTau.*fftGen(par.parphi.parFft,fft_u[par.order]))),0,par.parphi)))")

    result = zeros(typeof(par.parphi.eps), par.n_max+1, size(u0,1))
    result[1,:] = u0

    println("diff result ori = $(norm(u0-_getresult(fft_u[par.order], 0, par.parphi)))")

    _initAB(par, fftfct, fft_u)

    memfft = view(fftfct, par.order:(2par.order-1))

    for i=2:par.order
        result[i,:] = _getresult(fft_u[par.order+i-1], (i-1)*par.dt, par.parphi)
    end

     # ring permutation where the beginning becomes the end and the rest is shifted by one
    permut = collect(Iterators.flatten((2:par.order,1:1)))

    ut0_fft = fft_u[end]
    println("")
    for i=par.order:par.n_max
        if i%10000 == 1
            println(" $(i-1)/$(par.n_max)")
        end
        resfft, ut0_fft = _trAB(par, memfft, ut0_fft)
        result[i+1, :] = _getresult(ut0_fft, i*par.dt, par.parphi)
        memfft = memfft[permut]
        memfft[end] = resfft
        if i%100 == 0
            print("x")
        end
    end

    return result
end
