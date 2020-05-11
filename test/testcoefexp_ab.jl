include("../src/coefexp_ab.jl")
include("datacoefexp.jl")
include("datacoefexp10.jl")
using LinearAlgebra
using Random
using Test
function testcoefexp_ab()
    prec = prec_old = precision(BigFloat)
    @time @testset "test 1 coef for exponential Adams-Bashforth" begin
        for i=1:2
            setprecision(prec) do
                epsilon = parse(BigFloat,"0.1")
                dt = parse(BigFloat,"0.0001")
                par = CoefExpAB(15, epsilon, 32, dt)
                tab, list_j = get_coef_ab_for_test()
                for j in list_j
                    @time @test all(isapprox.(
    par.tab_coef[:,:,j], 
    tab[:,:,j], 
    atol=eps(BigFloat)*100, 
    rtol=eps(BigFloat)*100
))
                end
                @test par.tab_coef == -conj(par.tab_coef_neg)
            end
            prec *= 2
        end
    end
    setprecision(prec_old)
end
function testcoefexp10_ab()
    prec = prec_old = precision(BigFloat)
    @time @testset "test 1 coef for exponential Adams-Bashforth" begin
        for i=1:2
            setprecision(prec) do
                epsilon = big"1"/10^10
                dt = big"1"/10000
                par = CoefExpAB(15, epsilon, 32, dt)
                tab, list_j = get_coef_ab_for_test10()
                for j in list_j
                    @time @test all(isapprox.(
    par.tab_coef[:,:,j], 
    tab[:,:,j], 
    atol=eps(BigFloat)*100, 
    rtol=eps(BigFloat)*100
))
                end
                @test par.tab_coef == -conj(par.tab_coef_neg)
            end
            prec *= 2
        end
    end
    setprecision(prec_old)
end
function testpolylagrange()
     @time @testset "test polylagrange" begin
        polyOri = Polynomial(convert( Array{Complex{Rational{BigInt}}},
                                [3//4, 13//24, 1//7, 2//9, 40//51,
                                    1//128, 2//79, 3//334, 7//19, 1//100,
                                    1//144, 3//32, 7//19, 2//5, 1//101]
                                ) )
        s = size(polyOri,1)
        polyResult = Polynomial([zero(Complex{Rational{BigInt}})])
        for i = 1:s
            x = 1-i
            polyResult += polyOri(x)*getpolylagrange(i-1,s-1,BigInt)
        end
        @test polyOri == polyResult
    end
end
function testinterpolate()
    Random.seed!(9817171)
    @time @testset "testinterpolate" begin
        order= 8
        tab = Vector{Array{BigFloat,2}}(undef,order+1)
        for i=1:(order+1)
            tab[i]=zeros(BigFloat,4,32)
        end
        fctref = Array{Polynomial{BigFloat},2}(undef,4,32)
        for i=1:4 
            for j=1:32
                fctref[i,j]=Polynomial(rand(BigFloat,8)/10)
                for k = 1:(order+1) 
                    tab[k][i,j] = fctref[i,j](k-1)
                end
            end
        end
        value = big"4.22856371432981357654"
        res = interpolate(tab, order, value)
        dt = big"0.00123"
        t_deb=big"1.82898988989"
        value2=t_deb+value*dt
        tab_time=collect(t_deb:dt:(t_deb+(order+2)*dt))
        res2 = interpolate(tab_time, tab, order, value2)
        resref = zeros(BigFloat,4,32)
        for i=1:4 
            for j=1:32
                resref[i,j]=fctref[i,j](value)
            end
        end
        println("norm=$(norm(resref-res,Inf))")
        println("norm2=$(norm(resref-res2,Inf))")
        @test isapprox( resref, res,atol=1e-60,rtol=1e-60)
        @test isapprox( resref, res2,atol=1e-60,rtol=1e-48)
    end
end
testinterpolate()
testpolylagrange()
testcoefexp_ab()
testcoefexp10_ab()
