include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Test

function testtwoscales_pure_ab()

    u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]
    epsilon = big"0.000001"
    t_max = big"1.0"
    println("launching of julia solver")
    prob = ODEProblem(henon_heiles_julia, u0, (big"0.0", t_max), epsilon)
    @time sol = solve(prob , abstol=1e-17, reltol=1e-17)
    sol_ref = sol.u[end]
    nb=10
    order=3
    parphi = PreparePhi( 
        epsilon, 
        64, 
        [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], 
        henon_heiles
    )
    par_u0 = PrepareU0(parphi, order+2, u0)
    res_err = zeros(BigFloat,5)
    ind = 1
    while nb <= 100000
        pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
        result = twoscales_pure_ab(pargen)
        res_err[ind] = norm(result[:,end]-sol_ref,Inf)
        println("\nnb=$nb err=$(res_err[ind])")
        nb *= 10
        ind += 1
    end
end
testtwoscales_pure_ab()
