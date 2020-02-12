include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Test

function testtwoscales_pure_ab()

    u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]
    epsilon = big"0.01"
    t_max = big"1.0"
    prob = ODEProblem(henon_heiles_jul, u0, (big"0.0", t_max), epsilon)
    sol = solve(prob , abstol=1e-25, reltol=1e-25)
    sol_ref = sol.u[end]
    nb=10
    order=3

    res_err = zeros(BigFloat,3)
    ind = 1
    while nb <= 1000
        parphi = PreparePhi( 
    epsilon, 
    32, 
    [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], 
    henon_heiles
)
        par_u0 = PrepareU0(parphi, order+2, u0)   
        pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
        result = twoscales_pure_ab(pargen)
        res_err[ind] = norm(result[:,end]-sol_ref,Inf)
        println("nb=$nb err=$(res_err[ind])")
        nb *= 10
        ind += 1
    end
end
testtwoscales_pure_ab()
