include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")
using DifferentialEquations
using LinearAlgebra
using Random
using Test

function testtwoscales_pure_ab()

    u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]

    epsilon = big"0.01"
    t_max = big"1.0"
    println("launching of julia solver")
    prob = ODEProblem(henon_heiles_julia, u0, (big"0.0", t_max), epsilon)
    @time sol = solve(prob , abstol=1e-25, reltol=1e-25)
    sol_ref = sol.u[end]
    nb=10
    parphi = PreparePhi( 
        epsilon, 
        128, 
        [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], 
        henon_heiles
    )
    for order=2:10
        nb = 100
        par_u0 = PrepareU0(parphi, order+2, u0)
        res_err = zeros(BigFloat,5)
        ind = 1
        while nb <= 1000
            pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
            result = twoscales_pure_ab(pargen)
            res_err[ind] = norm(result[:,end]-sol_ref,Inf)
            println("\nnb=$nb order=$order err=$(res_err[ind])")
            nb *= 10
            ind += 1
        end
        coef = Float64(log(res_err[1]) - log(res_err[2]))/log(10)
        println("order=$order coef log = $coef")
    end
end




# function testtwoscales_pure_ab2()

#     u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]

#     B = [ big"0.12984599677" big"0.9277" big"0.32984110099677" big"0.142984599677"
#     big"0.4294599677" big"0.127337" big"0.4298411009977" big"0.99484599677"
#     big"0.2298499677" big"0.327667" big"0.1298410099677" big"0.342984599677"
#     big"0.7298459677" big"0.027887" big"0.7294110099677" big"0.66294599677"
#     ]
#     A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]

#     epsilon = big"0.1"

#     C = 1/epsilon * A + B


#     t_max = big"1.0"
#     fct_jul = (u,p,t) -> C*u
#     println("launching of julia solver")
#     prob = ODEProblem(fct_jul, u0, (big"0.0", t_max), epsilon)
#     @time sol = solve(prob , abstol=1e-20, reltol=1e-20)
#     sol_ref = sol.u[end]

#     sol_cal = exp(t_max*C)*u0
#     println("norm = $(norm(sol_ref-sol_cal))")
# end
function testtwoscales_pure_ab3()

    u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]

    Random.seed!(9988)

    # B = [ big"-0.12984599677" big"-0.9277" big"0.32984110099677" big"0.142984599677"
    # big"-0.4294599677" big"0.127337" big"0.4298411009977" big"0.99484599677"
    # big"0.2298499677" big"0.327667" big"0.1298410099677" big"-0.342984599677"
    # big"0.7298459677" big"-0.027887" big"0.7294110099677" big"-0.66294599677"
    # ]

    B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)

    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]

    epsilon = big"0.000001"

    C = 1/epsilon * A + B

    t_max = big"1.0"

    sol_ref = exp(t_max*C)*u0

    fct = u -> B*u

    parphi = PreparePhi( 
        epsilon, 
        32, 
        [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0], 
        fct,
        B
    )
    for order=2:10
        nb = 1000
        par_u0 = PrepareU0(parphi, order+2, u0)
        res_err = zeros(BigFloat,5)
        ind = 1
        while nb <= 10000
            pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
            result = twoscales_pure_ab(pargen)
            res_err[ind] = norm(result[:,end]-sol_ref,Inf)
            println("\nnb=$nb order=$order err=$(res_err[ind])")
            nb *= 10
            ind += 1
        end
        coef = Float64(log(res_err[1]) - log(res_err[2]))/log(10)
        println("order=$order coef log = $coef")
    end
end
function testtwoscales_pure_ab_epsilon()
    Random.seed!(19988)
    u0=rand(BigFloat,4)
    println("u0=$u0")
    B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    fct = u -> B*u
    order = 5
    for i=1:5
        epsilon = big"1.0"/big"10"^(i*3)
        eps_v = convert(Float32, epsilon)
        t_max = big"1.0"
        parphi = PreparePhi( 
            epsilon, 
            32, 
            A, 
            fct,
            B
        )    
        sol_ref = getexactsol(parphi, u0, t_max)
        println("epsilon=$eps_v sol_ref=$sol_ref")
        nb = 1000
        par_u0 = PrepareU0(parphi, order+2, u0)
        res_err = zeros(BigFloat,5)
        ind = 1
        while nb <= 10000
            pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
            result = twoscales_pure_ab(pargen)
            res_err[ind] = norm(result[:,end]-sol_ref,Inf)
            println("\nnb=$nb epsilon=$eps_v err=$(res_err[ind])")
            nb *= 10
            ind += 1
        end
        coef = Float64(log(res_err[1]) - log(res_err[2]))/log(10)
        eps_v
        println("epsilon=$eps_v coef log = $coef")
    end
end
function testtwoscales_interpolate()
    seed=9818871
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    order = 6
    nb = 100
    t_max = big"1.0"
    for i=1:10
        seed += 10
        Random.seed!(seed)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        fct = u -> B*u
        epsilon = big"0.4"/2.0^i
        eps_v = convert(Float32, epsilon)
        t_max = big"1.0"
        parphi = PreparePhi(epsilon, 32, A, fct, B)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
        @time result, tfft, tabu = twoscales_pure_ab(pargen, 
    only_end=false, res_fft=true)
        reftol=norm(getexactsol(parphi, u0, t_max)-result[:,end], Inf)*10
        for j=1:100
            t=rand(BigFloat)
            res_ex=getexactsol(parphi, u0, t)
            res_ap=_getresult(tabu, t, parphi, 0, t_max, order)
            println("t=$t norm=$(norm(res_ex-res_ap, Inf))")
            @test isapprox(res_ex, res_ap, atol=reftol, rtol=reftol*10)
        end
    end
end
function testtwoscales_short()
    seed=981885
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    order = 12
    nb = 10
    t_max = big"1.0"
    for i=1:10
        seed += 10
        Random.seed!(seed)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        fct = u -> B*u
        epsilon = big"0.00004"/2.0^i
        eps_v = convert(Float32, epsilon)
        t_max = big"0.1"
        parphi = PreparePhi(epsilon, 32, A, fct, B)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalePureAB(nb, t_max, order, par_u0)
        @time sol = twoscales_pure_ab(pargen, only_end=true)
        resnorm=norm(getexactsol(parphi, u0, t_max)-sol, Inf)
        println("resnorm=$resnorm")
        @test resnorm < 1e-10
    end
end
testtwoscales_short()
testtwoscales_interpolate()
testtwoscales_pure_ab()
# testtwoscales_pure_ab2()
testtwoscales_pure_ab3()
testtwoscales_pure_ab_epsilon()
