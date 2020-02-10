include("../src/preparephi.jl")
include("../src/henon_heiles.jl")
include("dataprep_u0.jl")
using Test
function testpreparephi()
    @time @testset "phi and prepare u0" begin

        epsilon = big"1"/256
        n_tau = 32
        parphi = PreparePhi(
    epsilon, 
    n_tau, 
    [0 0 1 0; 0 0 0 0; -1 0 0 0; 0 0 0 0],
    henon_heiles
)
        tab_ref = undef
        setprecision(1024) do
            tab_ref = get_prepare_u0_for_test()
        end
        u0 = BigFloat.([0.125, 0.140625, 0.15625, 0.171875])
        for ord=2:10
            paru0 = PrepareU0(parphi, ord, u0, 1024)
            println("norm=$(norm(tab_ref[:, :, ord]- paru0.ut0))")
            @test isapprox(tab_ref[:, :, ord], paru0.ut0, atol=1e-77, rtol=1e-77)
        end
    end
end
testpreparephi()
