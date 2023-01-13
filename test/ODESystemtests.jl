using DrWatson, Test
@quickactivate "tumor_control"
using Revise
using tumor_control: ODESystems

@testset "ODEsystems construction" begin

    p = (ρ=.035,
        m=30,
        K=4.8e6,
        α=0.01,
        β=0.1,
        C_max=5)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)

end

@testset "ODEsystems numeric integration" begin
    p = (ρ=.035,
        m=30,
        K=4.8e6,
        α=0.01,
        β=0.1,
        C_t=ODESystems.nullfunc)

    n_steps = 10
    Δts = repeat([.1],10)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)
    sol = ODESystems.solveEuler(system,[5e4,10],Δts) 
    println(sol)
    @test 1 == 1
end

@testset "ODEsystems JuMP control solving" begin
    #@test_broken
end