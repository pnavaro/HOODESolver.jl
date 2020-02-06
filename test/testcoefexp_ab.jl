include("../src/coefexp_ab.jl")
include("datacoefexp.jl")
using LinearAlgebra
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
function testpolylagrange()
     @time @testset "test polylagrange" begin
        polyOri = Poly(convert( Array{Complex{Rational{BigInt}}},
                                [3//4, 13//24, 1//7, 2//9, 40//51,
                                    1//128, 2//79, 3//334, 7//19, 1//100,
                                    1//144, 3//32, 7//19, 2//5, 1//101]
                                ) )
        s = size(polyOri,1)
        polyResult = Poly([zero(Complex{Rational{BigInt}})])
        for i = 1:s
            x = 1-i
            polyResult += polyOri(x)*getpolylagrange(i-1,s-1,BigInt)
        end
        @test polyOri == polyResult
    end
end
testpolylagrange()
testcoefexp_ab()
