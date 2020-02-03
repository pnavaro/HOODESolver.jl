
include("../src/expmatrices.jl")
using LinearAlgebra
using Random
using Test


function testexp()
    A=zeros(Float64, 4, 4)
    A[1,3] = 1
    A[3,1] = -1
    @time @testset "test exp big matrices" begin
        @test isapprox( one(A), _expm1(0A), atol=1e-16)
        @test isapprox(exp(0.123*A), _expm1(0.123*A), atol=1e-15)

        Abig = convert(Array{BigFloat,2},A)

        Random.seed!(81725)
        sens=1
        for i = 1:100:100000
            v=rand(BigFloat)*i*sens
            res = _expm1(v*Abig)
            resRef = one(Abig)
            resRef[1,1] = resRef[3,3] = cos(v)
            resRef[1,3] = sin(v)
            resRef[3,1]= -resRef[1,3]
            @test isapprox( resRef, res, atol=1e-76*i)
            sens *= -1
        end
        b = [ 0 0 0 1; 0 0 -1 0; 0 1 0 0;-1 0 0 0]
        for i=1:20
            C=rand(4,4)*i
            B = b*i
            @test isapprox(exp(C), _expm1(C), rtol=1e-14)
            @test isapprox(exp(10C), _expm1(10C), rtol=1e-13)
            @test isapprox(exp(100C), _expm1(100C), rtol=1e-12)
            @test isapprox(exp(1000C), _expm1(1000C), rtol=1e-11)
            @test isapprox(exp(10000C), _expm1(10000C), rtol=1e-10)
            @test isapprox(exp(B), _expm1(B), rtol=1e-14)
            @test isapprox(exp(10B), _expm1(10B), rtol=1e-13)
            @test isapprox(exp(100B), _expm1(100B), rtol=1e-12)
            @test isapprox(exp(1000B), _expm1(1000B), rtol=1e-11)
            @test isapprox(exp(10000B), _expm1(10000B), rtol=1e-10)

        end
    end
end

testexp()
