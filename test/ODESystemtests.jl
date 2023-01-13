using DrWatson, Test
@quickactivate "tumor_control"
using Revise
using tumor_control: ODESystems

@testset "ODEsystem construction" begin

    p = (ρ=.035,
        m=30,
        K=4.8e6,
        α=0.01,
        β=0.1,
        C_max=5)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)

end

@testset "ODEsystem numeric solving" begin

    @test 1 == 1
end