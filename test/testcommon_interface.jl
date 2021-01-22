using HOODESolver
using LinearAlgebra
using DoubleFloats
using Random
using Test
using DiffEqBase


function test_operator(T::DataType)
#    mat = [0 1;-1 0]
    for mat in [
    T.([0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]),
    T.([0 1;-1 0]),
    T.([0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0; 
     0 0 0 0 1 0; 0 0 0 -1 0 0; 0 0 0 0 0 0])
]
        epsilon = T(1//1024)
        linop = LinearHOODEOperator((1/epsilon)*mat)
        @test linop.epsilon == epsilon
        @test linop.A == mat
        linop2 = LinearHOODEOperator(linop)
        @test linop2.epsilon == epsilon
        @test linop2.A == mat
        linop3 = LinearHOODEOperator(epsilon, mat)
        @test linop3.epsilon == epsilon
        @test linop3.A == mat

    end
end


@testset "test operator" begin
    test_operator(Float64)
    test_operator(Double64)
    test_operator(BigFloat)
end

function fct4bis( du, u, p, t)
    B=[0.1 1.1 0.2 05 ; -0.9 -0.12 -0.7 0.4 ; 0.5 0.66 0.7 0.8 ; -0.34 0.8 0 0.3]
    du[:] = B*u
    missing
end

function testcommon_interface_fct(linop)
    B=[0.1 1.1 0.2 05 ; -0.9 -0.12 -0.7 0.4 ; 0.5 0.66 0.7 0.8 ; -0.34 0.8 0 0.3]
    u0 = 2rand(4)-ones(4)
    sol_ref = exp(1.0*(1/linop.epsilon*linop.A+B))*u0
    fct = (u,p,t) -> B*u
    prob = SplitODEProblem{false}(linop, fct, u0, (0.0, 1.0))
    sol = solve(prob, HOODEAB())
    @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
    fct = (u,p) -> B*u
    prob = SplitODEProblem{false}(linop, fct, u0, (0.0, 1.0))
    sol = solve(prob, HOODEAB())
    @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
    fct = (u) -> B*u
    prob = SplitODEProblem{false}(linop, fct, u0, (0.0, 1.0))
    sol = solve(prob, HOODEAB())
    @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
    # prob = SplitODEProblem(linop, fct4bis, u0, (0.0, 1.0), missing, A, epsilon)
    # sol = solve(prob)
    # @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
    
end

@time @testset "test interface type of function" begin
    A =  [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]
    epsilon = 0.0001
    linop = LinearHOODEOperator(epsilon,A)
    testcommon_interface_fct(linop)
end




