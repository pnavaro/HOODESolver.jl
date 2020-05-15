using LinearAlgebra
using Random
using Test
include("../src/interface.jl")
include("../src/henon_heiles.jl")


function fct4( du, u, p, t)
    B=[0.1 1.1 0.2 05 ; -0.9 -0.12 -0.7 0.4 ; 0.5 0.66 0.7 0.8 ; -0.34 0.8 0 0.3]
    du[:] = B*u
    missing
end
   

function testinterface_fct()
    @time @testset "test interface type of function" begin
        A =  [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0]
        B=[0.1 1.1 0.2 05 ; -0.9 -0.12 -0.7 0.4 ; 0.5 0.66 0.7 0.8 ; -0.34 0.8 0 0.3]
        u0 = 2rand(4)-ones(4)
        epsilon = 0.0001
        sol_ref = exp(1.0*(1/epsilon*A+B))*u0
        fct = (u,p,t) -> B*u
        prob = HiOscODEProblem(fct, u0, (0.0, 1.0), missing, A, epsilon)
        sol = solve(prob)
        @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
        fct = (u,p) -> B*u
        prob = HiOscODEProblem(fct, u0, (0.0, 1.0), missing, A, epsilon)
        sol = solve(prob)
        @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
        fct = (u) -> B*u
        prob = HiOscODEProblem(fct, u0, (0.0, 1.0), missing, A, epsilon)
        sol = solve(prob)
        @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
        prob = HiOscODEProblem(fct4, u0, (0.0, 1.0), missing, A, epsilon)
        sol = solve(prob)
        @test isapprox(sol_ref, sol[end], atol=1e-7, rtol=1e-6)
    end
end


function testinterface_epsilon()
    @time @testset "test interface while epsilon varying" begin
        Random.seed!(199881)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        fct = (u,p,t) -> B*u
        order = 5
        for i=1:1
            epsilon = big"1.0"/big"10"^(i*3)
            eps_v = convert(Float32, epsilon)
            t_max = big"1.0"
            sol_ref = exp(t_max*(1/epsilon*A+B))*u0
            prob = HiOscODEProblem(fct, u0, (big"0.0", t_max), missing, A, epsilon)
            println("epsilon=$eps_v sol_ref=$sol_ref")
            nb = 100
            res_err = zeros(BigFloat,5)
            ind = 1
            while nb <= 1000
                sol = solve(prob, nb_t=nb, order=5, dense=false)
                res_err[ind] = norm(sol[end]-sol_ref,Inf)
 #               println("\nnb=$nb epsilon=$eps_v err=$(res_err[ind])")
                nb *= 10
                ind += 1
            end
            coef = Float64(log(res_err[1]) - log(res_err[2]))/log(10)
            @test isapprox( order, coef, atol=0.01)
            println("epsilon=$eps_v coef log = $coef")
        end
    end
end
function testinterface_interpolate()
    @time @testset "test interface for interpolation" begin
        seed=123871
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        order = 6
        nb = 100
        t_max = big"1.0"
        for i=1:1
            seed = seed*1341+6576123
            Random.seed!(seed)
            u0=rand(BigFloat,4)
            println("u0=$u0")
            B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
            fct = (u,p,t) -> B*u
            epsilon = big"0.4"/2.0^i
            eps_v = convert(Float32, epsilon)
            t_max = big"1.0"
            prob = HiOscODEProblem(fct, u0, (big"0.0",t_max), missing, A, epsilon)
   #         sol = solve(prob, getprecision=false)
            sol = solve(prob)
            m = 1/epsilon*A+B
            reftol=norm(exp(t_max*m)*u0-sol[end], Inf)*10
            for j=1:10
                t=rand(BigFloat)
                res_ex=exp(t*m)*u0
                res_ap=sol(t)
                println("t=$t norm=$(norm(res_ex-res_ap, Inf))")
                @test isapprox(res_ex, res_ap, atol=reftol, rtol=reftol*10)
            end
        end
    end
end
function testinterface_interpolate_float()
    @time @testset "test interface for interpolation" begin
        seed=124333
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        order = 6
        nb = 100
        t_max = 1.0
        for i=1:1
            seed = seed*1341+6576123
            Random.seed!(seed)
            u0=rand(4)
            println("u0=$u0")
            B = 2rand(4,4)-ones(4,4)
            fct = (u,p,t) -> B*u
            epsilon = 0.4/2.0^i
            eps_v = convert(Float32, epsilon)
            t_max = 1.0
            prob = HiOscODEProblem(fct, u0, (0.0,t_max), missing, A, epsilon)
   #         sol = solve(prob, getprecision=false)
            sol = solve(prob)
            m = 1/epsilon*A+B
            reftol=norm(exp(t_max*m)*u0-sol[end], Inf)*10
            @test reftol < 1e-6
            for j=1:10
                t=rand()
                res_ex=exp(t*m)*u0
                res_ap=sol(t)
#                println("type ex=$(typeof(res_ex)) type ap=$(typeof(res_ap))")
#                println("t=$t norm=$(norm(res_ex-res_ap, Inf))")
                @test isapprox(res_ex, res_ap, atol=reftol, rtol=reftol*10)
            end
        end
    end
end
function testinterface_short()
    seed=981885
    A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
    order = 7
    nb = 5
     for i=1:10
        seed += 10
        Random.seed!(seed)
        u0=rand(BigFloat,4)
        println("u0=$u0")
        B = 2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        fct = (u,p,t) -> B*u
        epsilon = big"0.00004"/2.0^i
        eps_v = convert(Float32, epsilon)
        t_max = big"0.01"
        prob = HiOscODEProblem(fct, u0, (big"0.0",t_max), missing, A, epsilon)
        sol = solve(prob, getprecision=false, order=order, nb_t=nb)
        resnorm=norm(exp(t_max*(1/epsilon*A+B))*u0-sol[end], Inf)
        println("sol=$(sol[end])")
        println("solexact=$(exp(t_max*(1/epsilon*A+B))*u0)")
        println("resnorm=$resnorm")
        @test resnorm < 1e-10
    end
end
function tts_time(t_begin, t_end)
    @time @testset "test interface time from $t_begin to $t_end" begin
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        B =  2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        u0 = rand(BigFloat, 4)
        fct = (u,p,t) -> B*u
        epsilon=big"0.000000345"
        prob = HiOscODEProblem(fct, u0, (t_begin,t_end), missing, A, epsilon)
        sol = solve(prob, getprecision=false, nb_t=100, order=4)
        m = 1/epsilon*A+B
        solref = exp((t_end-t_begin)*m)*u0
        println("sol=$(sol[end]) solref=$solref norm=$(norm(sol[end]-solref,Inf))")
        @test isapprox(sol[end], solref,atol=1e-8, rtol=1e-7)
        for i=1:10
            t = rand(BigFloat)*(t_end-t_begin) + t_begin
            res_ex=exp((t-t_begin)*m)*u0
            res_ap=sol(t)
            @test isapprox(res_ex, res_ap, atol=1e-8, rtol=1e-7)
        end
    end
end
function tts_time_time(t_begin, t_end)
    @time @testset "test interface time time from $t_begin to $t_end" begin
        A =  [0 0 1 0; 0 0 0 0;-1 0 0 0; 0 0 0 0]
        B =  2rand(BigFloat,4,4)-ones(BigFloat,4,4)
        u0 = 2rand(BigFloat, 4)-ones(BigFloat,4)
        fct = (u,p,t) -> B*u +t*p[1] +p[2]
        tuple_p = (2rand(BigFloat,4)-ones(BigFloat,4),
    2rand(BigFloat,4)-ones(BigFloat,4))
        epsilon=big"0.000003467"
        nb = 100
        order = 4
        prob = HiOscODEProblem(fct, u0, (t_begin,t_end), tuple_p, A, epsilon, B)
        sol = solve(prob, getprecision=false, nb_t=100, order=4)
        solref = getexactsol(sol.par_u0.parphi, u0, t_end)
        println("sol=$(sol[end]) solref=$solref norm=$(norm(sol[end]-solref,Inf))")
        @test isapprox(sol[end], solref,atol=1e-8, rtol=1e-7)
        for i=1:10
            t = rand(BigFloat)*(t_end-t_begin) + t_begin
            res_ex=getexactsol(sol.par_u0.parphi, u0, t)
            res_ap=sol(t)
            @test isapprox(res_ex, res_ap, atol=1e-8, rtol=1e-7)
        end
    end
end
function testinterface_time()
    seed=7887
    Random.seed!(seed)
    tts_time(big"0.0",big"1.0")
    tts_time(big"0.0",big"1.4565656565")
    tts_time(big"2.234566",big"3.45766790123")
    tts_time(-big"0.8457676",big"0.56716")
    tts_time(-big"1.8111457676",-big"0.345456716")
    tts_time(big"1.8457676",big"0.56716")
    tts_time(-big"0.8457676",-big"1.56716")
    tts_time(-big"0.845722676",-big"0.56716")
end
function testinterface_time_time()
    seed=788227
    Random.seed!(seed)
    tts_time_time(big"0.0",big"1.0")
    tts_time_time(big"0.0",big"1.4565656565")
    tts_time_time(big"11.234566",big"12.45766790123")
    tts_time_time(-big"0.8457676",big"0.56716")
    tts_time_time(-big"1.8111457676",-big"0.345456716")
    tts_time_time(big"1.8457676",big"0.56716")
    tts_time_time(-big"0.8457676",-big"1.56716")
    tts_time_time(-big"0.845722676",-big"0.56716")
end

testinterface_fct()
testinterface_interpolate_float()
testinterface_time()
testinterface_time_time()
testinterface_epsilon()
testinterface_interpolate()
testinterface_short()
# testtwoscales_pure_ab()
# testtwoscales_pure_ab2()
#testtwoscales_pure_ab3()
