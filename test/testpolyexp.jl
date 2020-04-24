
include("../src/polyexp.jl")
using Test
function testPolyExp()

    @time @testset "test PolyExp" begin
        pol = PolyExp(convert( Vector{Complex{Rational{BigInt}}}, [2//3, 5//23, 57//4, 1//5]), big(im//7), big(1//10+0im) )

        c_res = coeffs(pol)

 #       println("pol=$(string(pol))")

        @test [2//3, 5//23, 57//4, 1//5, im//7, 1//10 ] == c_res

        polder = derivative(pol)

        c_resder = coeffs(polder)

        @test [2im//21+5//23, 5im//(23*7)+ 57//2, 57im//(4*7)+3//5, 1im//35,  im//7, 1//10 ] == c_resder

        polintder = integrate(polder)
        c_resintder = coeffs(polintder)

        @test c_res == c_resintder

        p1 = PolyExp([1.0, 2, 3], 3.4, 1.1)
        p2 = PolyExp([0.5, -1, 1], 1.3, 0.245)

        c_resmul = coeffs(p1*p2)
        @test isapprox([0.5, 0.0, 0.5, -1.0, 3.0, 4.7, 1.345], c_resmul, atol=1e-15)

        p3 = PolyExp([0.2, 0, 3.0, 2.0], 1.3, 0.245)

        c_resadd = coeffs(p2 + p3)
        @test isapprox([0.7, -1.0, 4.0, 2.0, 1.3, 0.245], c_resadd, atol=1e-15)
    end
end

testPolyExp()
