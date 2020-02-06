"""
    PreparePhi(nTau::Integer, eps::AbstractFloat)

Immutable structure, to share calculations, needed for the phi function.
These data can be used elsewhere for example in twoscale function.

# Arguments :
- `nTau::Integer` : number of value around the unit disk, it must be a power of two.
- `eps::AbstractFloat` : epsilon of the system.

# Implementation :
- tau : list of value of t from zero to 1 - 1/nTau.
- sinTauList : list of values sin( 2pi t ) for t in tau list.
- cosTauList : list of values cos( 2pi t ) for t in tau list.
- coefTauList : power coefficients.
- fftIntegralList : integral coefficents.
- parFft : fft parameters.
- tabA : differents values of A^n, with n ∈ N.

"""
struct PreparePhi
    epsilon
    n_tau
    list_tau
    matrix_A
    tau_A
    coefTauList
    fftIntegralList
    parFft
    tabA
    dictPhi
    nbCallPhi

    function  PreparePhi(n_tau::Integer, epsilon::AbstractFloat, matrix_A)
        @assert prevpow(2,n_tau) == n_tau "$nTau is not a power of 2"
        T = typeof(epsilon)

        tau_A = 
        tau_A = exp.([0:n_tau-1]*T(pi)/n_tau.*matrix_A)
        # Sinus and cosinus values around the unit circle.
        sinTauList = sinpi.(2tau)
        cosTauList = cospi.(2tau)

        coefTauList = [collect(0:nTau / 2 - 1); collect(-nTau / 2:-1)]

        # Coefficients to integrate
        fftIntegralList = -1im ./ coefTauList
        fftIntegralList[1] = 0

        if T==BigFloat
            parFft = PrepareFftBig(nTau,eps)
        else
            parFft = missing
        end

        A = zeros(Int64, 4, 4)

        tabA = Array{Int64,2}[]
        push!(tabA,A+1I)
        A[1,3] = 1
        A[3,1] = -1
        varA = A
        push!(tabA,varA)
        varA *= A
        while( varA != A )
            push!(tabA,varA)
            varA *= A
        end
        dictPhi=Array{Dict{Array{T,1},Array{T,2}},1}(undef,20)
        for i=1:20
            dictPhi[i] = Dict{Array{T,1},Array{T,2}}()
        end
        nbCallPhi = zeros(Int64, 20, 2)
        return new(
    nTau, eps, eps_rat,  tau, sinTauList, cosTauList, coefTauList, fftIntegralList, parFft, tabA, dictPhi, nbCallPhi
)

    end
end
