using DifferentialEquations
using LinearAlgebra
using Random
using Test
include("../src/twoscales_pure_ab.jl")
include("../src/henon_heiles.jl")

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
            pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
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

    fct = (u,p,t) -> B*u

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
            pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
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
    fct = (u,p,t) -> B*u
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
        nb = 100
        par_u0 = PrepareU0(parphi, order+2, u0)
        res_err = zeros(BigFloat,5)
        ind = 1
        while nb <= 1000
            pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
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
    for i=1:1
        seed += 10
        Random.seed!(seed)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        fct = (u,p,t) -> B*u
        epsilon = big"0.4"/2.0^i
        eps_v = convert(Float32, epsilon)
        t_max = big"1.0"
        parphi = PreparePhi(epsilon, 32, A, fct, B)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
        @time result, tfft, tabu = twoscales_pure_ab(pargen, 
    only_end=false, res_fft=true)
        reftol=norm(getexactsol(parphi, u0, t_max)-result[:,end], Inf)*10
        for j=1:10
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
    order = 7
    nb = 5
     for i=1:1
        seed += 10
        Random.seed!(seed)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        fct = (u,p,t) -> B*u
        epsilon = big"0.00004"/2.0^i
        eps_v = convert(Float32, epsilon)
        t_max = big"0.01"
        parphi = PreparePhi(epsilon, 32, A, fct, B)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalesPureAB(nb, t_max, order, par_u0)
        @time sol = twoscales_pure_ab(pargen, only_end=true)
        resnorm=norm(getexactsol(parphi, u0, t_max)-sol, Inf)
        println("resnorm=$resnorm")
        @test resnorm < 1e-10
    end
end
function tts_time(t_begin, t_end)
    @time @testset "test twoscales from $t_begin to $t_end" begin
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        B =  2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        u0 = rand(BigFloat, 4)
        fct = (u,p,t) -> B*u
        epsilon=big"0.00345"
        nb = 100
        order = 4
        ordprep= order+2
        parphi = PreparePhi(epsilon, 32, A, fct, B, t_0=t_begin)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalesPureAB(nb, t_end, order, par_u0)
        sol = twoscales_pure_ab(pargen, only_end=true)
        solref = getexactsol(parphi, u0, t_end)
        println("sol=$sol solref=$solref norm=$(norm(sol-solref,Inf))")
        @test isapprox(sol, solref,atol=1e-7, rtol=1e-6)
        result, tfft, tabu = twoscales_pure_ab(pargen, only_end=false, res_fft=true)
        for i=1:10
            t = rand(BigFloat)*(t_end-t_begin) + t_begin
            res_ex=getexactsol(parphi, u0, t)
            res_ap=_getresult(tabu, t, parphi, t_begin, t_end, order)
            println("i=$i")
            println("t=$t")
            println("res_ex=$res_ex")
            println("res_ap=$res_ap")
            @test isapprox(res_ex, res_ap, atol=1e-6, rtol=1e-5)
        end
    end
end
function tts_time_time(t_begin, t_end)
    @time @testset "test twoscales time time from $t_begin to $t_end" begin
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        B =  2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        u0 = 2rand(BigFloat, 4)-ones(BigFloat,4)
        fct = (u,p,t) -> B*u +t*p[1] +p[2]
        tuple_p = (2rand(BigFloat,4)-ones(BigFloat,4),
    2rand(BigFloat,4)-ones(BigFloat,4))
        epsilon=big"0.0003467"
        nb = 100
        order = 4
        ordprep= order+2
        parphi = PreparePhi(epsilon, 32, A, fct, B, paramfct=tuple_p, t_0=t_begin)
        par_u0 = PrepareU0(parphi, order+2, u0)
        pargen = PrepareTwoScalesPureAB(nb, t_end, order, par_u0)
        sol = twoscales_pure_ab(pargen, only_end=true)
        solref = getexactsol(parphi, u0, t_end)
        println("sol=$sol solref=$solref norm=$(norm(sol-solref,Inf))")
        @test isapprox(sol, solref,atol=1e-8, rtol=1e-7)
        result, tfft, tabu = twoscales_pure_ab(pargen, only_end=false, res_fft=true)
        for i=1:10
            t = rand(BigFloat)*(t_end-t_begin) + t_begin
            res_ex=getexactsol(parphi, u0, t)
            res_ap=_getresult(tabu, t, parphi, t_begin, t_end, order)
            @test isapprox(res_ex, res_ap, atol=1e-7, rtol=1e-6)
        end
    end
end
function testtwoscales_time()
    seed=7887
    Random.seed!(seed)
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    order = 5
    nb = 10
    t_max = big"1.0"
    tts_time(big"0.0",big"1.0")
    tts_time(big"0.0",big"1.4565656565")
    tts_time(big"2.234566",big"3.45766790123")
    tts_time(-big"0.8457676",big"0.56716")
    tts_time(-big"1.8111457676",-big"0.345456716")
    tts_time(big"1.8457676",big"0.56716")
    tts_time(-big"0.8457676",-big"1.56716")
    tts_time(-big"0.845722676",-big"0.56716")
end
function testtwoscales_time_time()
    seed=788227
    Random.seed!(seed)
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    order = 5
    nb = 10
    t_max = big"1.0"
    tts_time_time(big"0.0",big"1.0")
    tts_time_time(big"0.0",big"1.4565656565")
    tts_time_time(big"11.234566",big"12.45766790123")
    tts_time_time(-big"0.8457676",big"0.56716")
    tts_time_time(-big"1.8111457676",-big"0.345456716")
    tts_time_time(big"1.8457676",big"0.56716")
    tts_time_time(-big"0.8457676",-big"1.56716")
    tts_time_time(-big"0.845722676",-big"0.56716")
end
function test_exact()
    A = [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]

    epsilon = big"0.1"
    n_tau = 32
    parphi = PreparePhi(
        epsilon, 
        n_tau, 
        A,
        henon_heiles
    )
    @test !isexactsol(parphi)

    u0=[big"0.12345678", big"0.13787878", big"0.120099999001", big"0.12715124"]

    B = [ big"-0.12984599677" big"-0.9277" big"0.32984110099677" big"0.142984599677"
    big"-0.4294599677" big"0.127337" big"0.4298411009977" big"0.99484599677"
    big"0.2298499677" big"0.327667" big"0.1298410099677" big"-0.342984599677"
    big"0.7298459677" big"-0.027887" big"0.7294110099677" big"-0.66294599677"
    ]

    fct = u -> B*u

    parphi = PreparePhi(
        epsilon, 
        n_tau, 
        A,
        fct,
        B
    )
    @test isexactsol(parphi)

    C = 1/epsilon * A + B

    t = big"0.12347171717"

    @test  isapprox(exp(t*C)*u0, getexactsol(parphi, u0, t), atol=1e-77,rtol=1e-77)

end   
test_exact()
testtwoscales_time()

testtwoscales_time_time()
testtwoscales_short()
testtwoscales_interpolate()
# testtwoscales_pure_ab()
# testtwoscales_pure_ab2()
#testtwoscales_pure_ab3()
testtwoscales_pure_ab_epsilon()
