using NLControl: ODESystems
using JuMP
using Ipopt

@testset "StepFunction construction and use" begin
    # Test inner constructor
    @test_throws "Mismatched" ODESystems.StepFunc1DData([0.0,1.0,2.0],[0.0,1.0,2.0])
    @test_throws "Only 1D" ODESystems.StepFunc1DData(zeros(2,2),zeros(3,2))
    @test_throws "sorted" ODESystems.StepFunc1DData([3.5,1.0],[0.0, 0.0, 0.0])
    @test_throws "unique" ODESystems.StepFunc1DData([1.0,1.0],[0.0, 0.0, 0.0])

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
            (:C_data,ODESystems.StepFunc1DData(0.0))])

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!,2)

end

@testset "ODEsystems numeric integration" begin

    p = Dict([(:ρ,.035),
            (:m,30),
            (:K,4.8e6),
            (:α,0.01),
            (:β,0.1),
            (:C_max,5),
            (:C_data,ODESystems.StepFunc1DData(0.0))])

    n_steps = 10
    Δts = repeat([.1],10)

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!, 2)
    sol = ODESystems.solveEuler(system,[5e4,10],Δts) 
    #println(sol)

    sol = ODESystems.solve_diffeq(system,[5e4,10],Δts) 
    #println(sol)
    #println(transpose(reduce(hcat,sol.u)))
end

@testset "ODEsystems JuMP control variable creation" begin
    using JuMP
    jumpMod = Model()
    control_data = ODESystems.createJuMPControlVariable!(jumpMod,
                                                         :C_t)
    println(control_data)
    show(jumpMod)
    # TODO make more rigorous testing with different situations
end

@testset "ODEsystems JuMP initialization and solving" begin

    jumpMod = Model(Ipopt.Optimizer)
    control_data = ODESystems.createJuMPControlVariable!(jumpMod,
                                                         :C_t,
                                                         upper_bound=5)
    p = Dict([(:ρ,.035),
    (:m,30),
    (:K,4.8e6),
    (:α,0.01),
    (:β,0.1),
    (:C_max,5),
    (:C_data,control_data)])

    system = ODESystems.DynamicalSystem(p, ODESystems.twoPopRSC!, 2)

    n_steps = 10
    Δts = repeat([.1],100)
    jumpMod,u,du,t = ODESystems.addJuMPNLDynamics!(jumpMod,system, Δts, [5e4,10],[:C_data])
    @objective(jumpMod,Min,sum(u[:,end].^2))
    #show(jumpMod)
    #println()
    #@profview optimize!(jumpMod)
    optimize!(jumpMod)
end
