"""
Collection of functions implementing ODE dynamics using DifferentialEquations.jl style.
"""
module ODESystems
    #export twoPopRSC!
    # check this package out: https://juliapackages.com/p/parameterizedfunctions
    using UnPack
    using DifferentialEquations
    using ..tumor_control: NumInt
    import JuMP
    using JuMP: Model,@variables,fix,@NLexpressions,@NLconstraint
    using Ipopt

    """
    Wrapper class to interface tools from JuMP and SciML.
    """
    struct DynamicalSystem
        p::NamedTuple
        # C_t::Function{t::Number}
        dynamics!::Function # f!(du,u,p,t) -> du
        # stateVarNames::Dict
    end

    # TODO Create a system solution struct?
    # enables creating plot methods, compare methods

    # TODO: check if this is actually creating performance gains
    function generateCopyDynamics(system::DynamicalSystem)::Function
        function dynamics(du,u,p,t)
            du = similar(u)
            system.dynamics(du,u,p,t)
            return du
        end
        return dynamics
    end

    function solve_diffeq(system::DynamicalSystem,
                            u0::AbstractArray{T},
                            Δt::AbstractArray{T};
                            kwargs...) where T<:Number
        tspan = (0.0,sum(Δt))
        savetimes = append!([0.0],cumsum(Δt))
        ode = DifferentialEquations.ODEProblem(system.dynamics!,u0,tspan,system.p)
        solution = DifferentialEquations.solve(ode; kwargs..., saveat=savetimes)
        return solution
    end

    function solveEuler(system::DynamicalSystem,
                        u0::AbstractArray{T},
                        Δt::AbstractArray{T}) where T<:Number        
        return NumInt.eulerSolve(system.dynamics!,system.p,u0,Δt)
    end

    function solveRK4()
    #TODO: implement RK4
    end

    """

    Create a base JuMP model based on the provided dynamical system.
    The returned model object only contains dynamics and starting solutions
    """
    function createJuMPNLControlModel(system::DynamicalSystem, Δt, u0,solver=Ipopt.Optimizer)
        # Using user defined functions:
        # check https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/tips_and_tricks/#User-defined-functions-with-vector-outputs
        # https://jump.dev/JuMP.jl/stable/manual/nlp/#User-defined-Functions

        # Instantiate a JuMP model
        model = Model(solver)
        # Compute initializing solution
        u_init = solveEuler(system,u0,Δt)
        # Initialize it 
        @variables(model,begin
            u[i=1:n,j=1:m] ≥ 0, (start=u_init[i,j])
        end)
        
        # Create a time variable
        fix(t[i=1:n], append!([0.0],cumsum(Δt))[i])
        # Fix initial state
        fix(u[0,:], u0; force = true)

        # Create an expression containing the dynamics 
        # TODO: check that this works using an RK4 integration
        copydyn = generateCopyDynamics(system)
        @NLexpressions(model,begin
            # du/dt
            du[j = 1:n], copydyn(u[j],system.p,t[j])
        end)

        # Add the dynamics constraint based on the integration scheme
        @NLconstraint(model, [j=2:n], u[j] == eulerStep(u[j - 1], du[j - 1], Δt[j-1]))

        return model
    end

    """
    The control term should be more coarse grained than the time step used for solving.
    To avoid stiffness issues the control times will be rounded to the closest numerical solver timesteps at solving time.
    """
    function addJuMPControlArray()

    end

    function solveJuMPNLControlModel()
        
    end

    function fitLS()
        
    end

    function fitstiffLS()
        
    end

    function fitBayes()
        
    end

    """
        twoPopRSC!(du,u,p,t)

    Implements an ODE system for two cell populations (sensitive and resistant) 
    and a control term varying across time.
    
    ```math
    \\begin{cases}
    \\frac{ds}{dt} = ρ.s(t).(1-\\frac{s(t)+m.r(t)}{K}) - αC(t)s(t) \\\\
    \\frac{dr}{dt} = ρ.r(t).(1-\\frac{s(t)+m.r(t)}{K}) - β \\frac{s(t)r(t)}{K} \\\\
    \\end{cases}
    ```

    **Arguments**

    p: a named tuple/Dict or Struct containing the following symbols:
    - ρ: scalar
    - m: scalar
    - K: scalar
    - α: scalar
    - β: scalar
    - C_t: a 1 argument function returning the control value at time ``t``.

    """
    function twoPopRSC!(du,u,p,t)
        s,r = u
        # Properly unpack parameters: https://stackoverflow.com/questions/44298860/julia-best-practice-to-unpack-parameters-inside-a-function
        @unpack ρ,m,K,α,β,C_t = p

        # Sensitive cells dynamics
        du[1] = s′ = ρ*s*(1-((s+m*r)/K)) - α*C_t(t)*s
        # Resistant cells dynamics
        du[2] = r′ = ρ*r*(1-((s+m*r)/K)) - β*r*s/K
    end

    """
        constfunc_generator(c::Number)
    Generator to build one dimensionnal constant scalar function.

    Return a function object returning ``c`` for any passed value ``t``
    """
    constfunc_generator(c::Number) = (t::Number->c)

    """
        null_func(t::Number)
    Dummy function for a null control over any time ``t``.

    Can be used to easily solve a system containing a control term in a null control situation.
    """
    nullfunc(t::Number) = 0.0


end