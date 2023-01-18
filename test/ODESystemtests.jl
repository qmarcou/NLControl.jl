using DrWatson, Test
@quickactivate "tumor_control"
using Revise
using tumor_control: ODESystems

@testset "StepFunction construction and use" begin
    # Test inner constructor
    @test_throws "Mismatched" ODESystems.StepFunc1DData([0.0,1.0,2.0],[0.0,1.0,2.0])
    @test_throws "Only 1D" ODESystems.StepFunc1DData(zeros(2,2),zeros(3,2))
    @test_throws "sorted" ODESystems.StepFunc1DData([3.5,1.0],[0.0, 0.0, 0.0])

    # Test basic functionnality
    stepdata = ODESystems.StepFunc1DData([1.0,2.0],[0.0,1.0,2.0])
    @test ODESystems.stepFuncVal(stepdata,-5) == 0.0
    @test ODESystems.stepFuncVal(stepdata,1.0) == 1.0
    @test ODESystems.stepFuncVal(stepdata,1.5) == 1.0
    @test ODESystems.stepFuncVal(stepdata,2.5) == 2.0

    # Test outer constructors
    stepdata = ODESystems.StepFunc1DData(5.0)
    @test ODESystems.stepFuncVal(stepdata,-5) == 5.0
    @test ODESystems.stepFuncVal(stepdata,5) == 5.0

end

@testset "ODEsystems construction" begin

    p = Dict([(:ρ,.035),
            (:m,30),
            (:K,4.8e6),
            (:α,0.01),
            (:β,0.1),
            (:C_max,5),
            (:C_t,ODESystems.nullfunc)])

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)

end

@testset "ODEsystems numeric integration" begin

    p = Dict([(:ρ,.035),
            (:m,30),
            (:K,4.8e6),
            (:α,0.01),
            (:β,0.1),
            (:C_max,5),
            (:C_t,ODESystems.nullfunc)])

    n_steps = 10
    Δts = repeat([.1],10)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)
    sol = ODESystems.solveEuler(system,[5e4,10],Δts) 
    #println(sol)

    sol = ODESystems.solve_diffeq(system,[5e4,10],Δts) 
    #println(sol)
    println(transpose(reduce(hcat,sol.u)))
end

@testset "ODEsystems JuMP control solving" begin
    p = (ρ=.035,
    m=30,
    K=4.8e6,
    α=0.01,
    β=0.1,
    C_t=ODESystems.nullfunc)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!)

    n_steps = 10
    Δts = repeat([.1],10)
    jumpMod = ODESystems.createJuMPNLControlModel(system, Δts, [5e4,10])
end