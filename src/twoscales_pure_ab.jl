
include("preparephi.jl")
include("coefexp_ab.jl")
"""
    PrepareTwoScalesPureAB(nb_t, t_max, order, par_u0::PrepareU0)

Immutable structure, to share calculations, needed for the twoscale function

# Arguments :
- `nb_t::Int` : number of time slices
- `t_max`: end of the time
- `order` : order for compute the coefficients
- `par_u0::PrepareU0` : prepared initial data

# Keywords :
- `p_coef::Union{CoefExpAB,Missing}=missing` : precomputed coefficients of AB method 
- `verbose=100` : trace level

# Fields :
- nb_t : number of time slices
- t_max : end of the time
- order : order for compute the coefficients
- parphi : prepared parameters for phi (from par_u0)
- par_u0 : prepared initial data
- p_coef : computed coefficients
- exptau : exp( -im*dt*'\tau'/epsilon) for all '\tau' values
- exptau_inv : inverse of exptau
- verbose : trace level

"""
struct PrepareTwoScalesPureAB
    nb_t
    t_max
    dt
    order
    parphi::PreparePhi
    par_u0::PrepareU0
    p_coef::CoefExpAB
    exptau
    exptau_inv
    verbose   
    function PrepareTwoScalesPureAB(nb_t, t_max, order, par_u0::PrepareU0;
    p_coef::Union{CoefExpAB,Missing}=missing,
    verbose=100
)
        parphi = par_u0.parphi
        T = typeof(parphi.epsilon)
        dt = T(t_max-parphi.t_0)/nb_t
        p_coef = if ismissing(p_coef)
            CoefExpAB(order, parphi.epsilon, parphi.n_tau, dt)
        else
            p_coef
        end
        exptau = collect(transpose(exp.(-im*dt / parphi.epsilon * parphi.tau_list)))
        exptau_inv = collect(transpose(exp.(im*dt / parphi.epsilon * parphi.tau_list)))
        return new(
    nb_t,
    t_max, 
    dt, 
    order, 
    parphi,
    par_u0,
    p_coef, 
    exptau,
    exptau_inv,
    verbose
)
    end
end

function _calculfft(parphi::PreparePhi, u_caret)
    f = filtredfct(parphi, real(ifftgen(parphi.par_fft, u_caret)))
    return fftgen(parphi.par_fft, f)
end

function _calcul_ab(par::PrepareTwoScalesPureAB, ord, fftfct, fft_u, dec, sens)
    resfft = fft_u[dec-sens] .* ((sens==1) ? par.exptau : par.exptau_inv)
    tab_coef = (sens == 1) ? par.p_coef.tab_coef : par.p_coef.tab_coef_neg
    for k=1:ord
        indice = dec-sens*k
        resfft .+= transpose(tab_coef[:, k, ord]).*fftfct[indice]
    end
    fft_u[dec] = resfft
    if par.verbose >= 200 && isexactsol(par.parphi)
        t = par.parphi.t_0 + (dec-par.order)*par.dt
        u_exact = getexactsol(par.parphi, par.par_u0.up0[1:(end-1)], t)
        u_computed = _getresult( fft_u[dec], t-par.parphi.t_0, par.parphi)
        err = Float64(norm(u_exact-u_computed, Inf))
        i=dec-par.order
        traceln(200, "i=$i err=$err", verbose=par.verbose)
    end
    fftfct[dec] = _calculfft(par.parphi, resfft)
end

# function _init_ab(par::PrepareTwoScalesPureAB, fftfct, fft_u)
# #    println("_init_ab order=$(par.order)")
#     if par.order != 1
#         for new_ord=2:par.order
#             _calcul_ab(par, new_ord-1, fftfct, fft_u, par.order+new_ord-1, 1)
#             for k = 1:new_ord-1
#                 _calcul_ab(par, new_ord, fftfct, fft_u, par.order-k, -1)
#             end
#             for k = 1:new_ord-1
#                 _calcul_ab(par, new_ord, fftfct, fft_u, par.order+k, 1)
#             end
#         end
#     end
# end
function _init_ab(par::PrepareTwoScalesPureAB, fftfct, fft_u)
    #    println("_init_ab order=$(par.order)")
        if par.order != 1
            for new_ord=2:par.order
                for k = 1:new_ord-1
                    _calcul_ab(par, new_ord-1, fftfct, fft_u, par.order-k, -1)
                end
                for k = 1:new_ord-1
                    _calcul_ab(par, new_ord, fftfct, fft_u, par.order+k, 1)
                end
            end
        end
    end
    
function _tr_ab(par::PrepareTwoScalesPureAB, fftfct, u_chap)
    resu_c = par.exptau.* u_chap
    bound = par.order-1
    for k =0:bound
        resu_c .+= transpose(par.p_coef.tab_coef[:, k+1, par.order]).*fftfct[end-k]
    end
    resfft = _calculfft(par.parphi, resu_c)
    return resfft, resu_c
end
function traceln( refniv, str; verbose=100)::Nothing
    if verbose >= refniv
        println(str)
    end
end
function trace( refniv, str; verbose=100)::Nothing
    if verbose >= refniv
        print(str)
    end
end

# for this function only, t is the time from the beginning
function _getresult(u_chap, t, par::PreparePhi)
    # matlab : u1=real(sum(fft(ut)/Ntau.*exp(1i*Ltau*T/ep)));
    u1 = real(u_chap * exp.(1im * par.tau_list * t / par.epsilon)) / par.n_tau
    res = exp( t / par.epsilon *par.sparse_Ap) * u1
    return res[1:(end-1)]
end
function _getresult( tab_u_chap, t, par::PreparePhi, t_begin, t_max, order)
    nb = size(tab_u_chap,1)-1
    dt = (t_max-t_begin)/nb
    t_ex = (t-t_begin)/dt
    t_int = convert(Int64,floor(t_ex))
    t_int_begin = t_int-div(order,2)
    if t_int_begin < 0
        t_int_begin = 0
    end
    if t_int_begin > nb-order
        t_int_begin = nb-order
    end
    t_ex -= t_int_begin
    t1 = t_int_begin+1
    N = (typeof(par.epsilon)==BigFloat) ? BigInt : Int64
    u_chap = interpolate(tab_u_chap[t1:(t1+order)], order, t_ex, N)
    return _getresult( u_chap, t-t_begin, par)
end
function _getresult( tab_t, tab_u_chap, t, par::PreparePhi, t_begin, t_max, order)
    nb = size(tab_u_chap,1)-1
    dt = (t_max-t_begin)/nb
    t_ex = (t-t_begin)/dt
    t_int = convert(Int64,floor(t_ex))
    t_int_begin = t_int-div(order,2)
    if t_int_begin < 0
        t_int_begin = 0
    end
    if t_int_begin > nb-order
        t_int_begin = nb-order
    end
    t1 = t_int_begin+1
    u_chap = interpolate(tab_t[t1:(t1+order)], tab_u_chap[t1:(t1+order)], order, t)
    return _getresult( u_chap, t-t_begin, par)
end
"""
    twoscales_pure_ab(par::PrepareTwoScalesPureAB; only_end::Bool=false, diff_fft::Bool=false, res_fft::Bool=false, verbose::Integer=100)

compute the data to get solution of the differential equation

# Arguments :
- `par::PrepareTwoScalesPureAB` : contains all the parameters and prepared data

# Keywords :
- `only_end=false` : if true return only the result for t_end
- `diff_fft::Bool=false` : if true return data about diff
- `res_fft::Bool=false` : if true return u_caret data indispensable for interpolation
- `verbose::Integer`: level off traces (0 means no output)

"""
function twoscales_pure_ab(
    par::PrepareTwoScalesPureAB; 
    only_end::Bool=false, diff_fft::Bool=false, 
    res_fft::Bool=false, verbose::Integer=100
)
    levelref=90
    T = typeof(par.parphi.epsilon)
    fftfct = Vector{Array{Complex{T}, 2}}(undef, 2par.order-1)
    fft_u = Vector{Array{Complex{T}, 2}}(undef, 2par.order-1)
    res_u = par.par_u0.ut0
    up0 = par.par_u0.up0
    fftfct[par.order] = fftgen(
    par.parphi.par_fft, 
    filtredfct(par.parphi, res_u)
    )
    fft_u[par.order] = fftgen(par.parphi.par_fft, res_u)
    result = zeros(typeof(par.parphi.epsilon), par.parphi.size_vect-1, par.nb_t+1)
    result[:,1] = up0[1:(end-1)]
    result_fft = undef
    _init_ab(par, fftfct, fft_u)
    memfft = fftfct[par.order:(2par.order-1)]
    if res_fft
        result_fft = Vector{Array{Complex{T}, 2}}(undef, par.nb_t+1)
        res_u_chap = Vector{Array{Complex{T}, 2}}(undef, par.nb_t+1)
        for i=1:min(par.order,par.nb_t+1)
            result_fft[i] = memfft[i]
        end
        for i=1:min(par.order,par.nb_t+1)
            res_u_chap[i] = fft_u[i+par.order-1]
        end
    end
    if !only_end
        for i=2:min(par.order,par.nb_t+1)
            result[:,i] = _getresult(fft_u[par.order+i-1], (i-1)*par.dt, par.parphi)
        end
    end
    tabdifffft = undef
    tabdifffft_2 = undef
    if diff_fft
        tabdifffft = zeros(BigFloat, par.nb_t)
        tabdifffft_2 = zeros(BigFloat, par.nb_t)
        for i = 1:(par.order-1)
            tabdifffft[i] = norm(memfft[i]-memfft[i+1],Inf)
            tabdifffft_2[i] = norm(memfft[i]-memfft[i+1])
        end
    end
     # ring permutation where the beginning becomes the end and the rest is shifted by one
    permut = collect(Iterators.flatten((2:par.order,1:1)))
    ut0_fft = fft_u[par.order-1+min(par.order,par.nb_t+1)]
    traceln(levelref, "", verbose=verbose)
    norm_delta_fft = 0
    norm_delta_fft_2 = 0
    nbnan = 0
    borne_nm=0
    c_mult=1.1
    for i=par.order:par.nb_t
        resfft, ut0_fft = _tr_ab(par, memfft, ut0_fft)
        if res_fft
            result_fft[i] = resfft
#            res_u_chap[i+par.order] = ut0_fft
            res_u_chap[i+1] = ut0_fft
        end
        if !only_end
            result[:,i+1] = _getresult(ut0_fft, i*par.dt, par.parphi)
        end
        if diff_fft
            nm = norm(resfft-memfft[end],Inf)
            nm_2 = norm(resfft-memfft[end])
            norm_delta_fft = max(nm,norm_delta_fft)
            norm_delta_fft_2 = max(nm_2,norm_delta_fft_2)
            tabdifffft[i] = nm
            tabdifffft_2[i] = nm_2
        end
        memfft = memfft[permut]
        memfft[end] = resfft
        if i%100 == 0
            trace(levelref, "x", verbose=verbose)
#            GC.gc()
        end
        if i%10000 == 0 || i == par.nb_t
            traceln(levelref, " $i/$(par.nb_t)", verbose=verbose)
        end
    end

#    println("norm diff fft = $norm_delta_fft")

    ret =  only_end ? _getresult(ut0_fft, par.t_max-par.parphi.t_0, par.parphi) : result

    if res_fft
        if diff_fft 
            return ret, result_fft, res_u_chap, tabdifffft, tabdifffft_2, norm_delta_fft, norm_delta_fft_2
        else
            return ret, result_fft, res_u_chap
        end
    else
        if diff_fft 
            return ret, tabdifffft, tabdifffft_2, norm_delta_fft, norm_delta_fft_2
        else
            return ret
        end
    end
end
