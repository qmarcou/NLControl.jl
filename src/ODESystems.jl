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
    using JuMP: Model,@variable,@variables,fix,@NLexpressions,@NLconstraint
    using Ipopt

    """
    An object storing information to build a 1D piecewise constant function.
    """
    struct StepFunc1DData
        cutpoints::AbstractArray
        values::AbstractArray
        #TODO add bounds on function space
        function StepFunc1DData(cutpoints,values)
            if ndims(values)!=1 | ndims(cutpoints)!=1
                error("Only 1D stepfunctions supported.")
            elseif length(cutpoints) != length(values)-1
                error("Mismatched size values and cutpoints arrays.")
            elseif (sort(cutpoints)!=cutpoints) | (unique(cutpoints)!=cutpoints)
                error("Cutpoints must be sorted and unique.")
            else
                new(cutpoints,values)
            end
        end
    end

    function StepFunc1DData(value)
        return StepFunc1DData(zeros(0),[value])
    end

#=     function StepFuncData(value,ncuts::Integer)
        
    end

    function StepFuncData()
        
    end =#

    function stepFuncVal(StepFunc1DData, x)
        index = findfirst(StepFunc1DData.cutpoints.>x)
        if isnothing(index)
            return last(StepFunc1DData.values)
        else
            return StepFunc1DData.values[index]
        end
    end

    """
    Wrapper class to interface tools from JuMP and SciML.
    """
    struct DynamicalSystem
        p::Dict{Symbol,Any}
        dynamics!::Function # f!(du,u,p,t) -> du
        ndims::Int64 # Number of u dimensions
        # stateVarNames::Dict
        function DynamicalSystem(p,dynamics!,ndims)
            if ndims<=0
                error("Dynamical system cannot be less than 1 dimensional.")
            else
                # TODO check that dynamics and ndims match
                return new(p,dynamics!,ndims)
            end
        end
    end

    # TODO: check if this is actually creating performance gains
    function generateCopyDynamics(system::DynamicalSystem)::Function
        function dynamics(du,u,p,t)
            du = similar(u)
            system.dynamics(du,u,p,t)
            return du
        end
        return dynamics
    end

    function generateSplattedScalarOutputCopyDynamics(system::DynamicalSystem)::Function
        nothing
    end

    function solve_diffeq(system::DynamicalSystem,
                            u0::AbstractArray{T},
                            Δt::AbstractArray{T};
                            kwargs...) where T<:Number
        tspan = (0.0,sum(Δt))
        savetimes = append!([0.0],cumsum(Δt))
        ode = DifferentialEquations.ODEProblem(system.dynamics!,u0,tspan,NamedTuple(system.p))
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

#=     """

    Create a base JuMP model based on the provided dynamical system.
    The returned model object only contains dynamics and starting solutions
    """
    function addJuMPNLDynamics(model::Model,
        system::DynamicalSystem,
        Δt,
        u0,
        controlVarNames::AbstractArray{Union{String,Symbol}})
        # Using user defined functions:
        # check https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/tips_and_tricks/#User-defined-functions-with-vector-outputs
        # https://jump.dev/JuMP.jl/stable/manual/nlp/#User-defined-Functions

        # Get array sizes
        n_rows=length(u0)
        n_cols = length(Δt) + 1 

        # Make a working copy of the parameters Dict
        p_copy = copy(system.p)

        for var in controlVarNames
            symbolic_val = p_copy[var]
            p_copy[var] = StepFunc1DData(symbolic_val.cutpoints,start_value.(symbolic_val.values))
        end

        # Make a copy of the system containing the parameters copy
        system_copy = DynamicalSystem(p_copy,system.dynamics!)

        # Compute initializing solution using the previously provided control starting values
        u_init = solveEuler(system_copy,u0,Δt)

        # Initialize it 
        @variables(model,begin
            u[i=1:n_rows,j=1:n_cols] ≥ 0, (start=u_init[i,j])
            t[i=1:n_cols] ≥ 0 # create a time vector variable 
        end)
        
        # Create a fixed time variable
        fix(t[i=1:n], append!([0.0],cumsum(Δt))[i])
        # Fix initial state to u0
        fix(u[0,:], u0; force = true)

        # Create an expression containing the dynamics 
        # TODO: check that this works using an RK4 integration
        copydyn = generateCopyDynamics(system)
        @NLexpressions(model,begin
            # du/dt
            du[j = 1:n_cols], copydyn(u[j],system.p,t[j])
        end)

        # Add the dynamics constraint based on the integration scheme
        @NLconstraint(model, [j=2:n], u[j] == eulerStep(u[j - 1], du[j - 1], Δt[j-1]))

        return model
    end =#

    """
    The control term should be more coarse grained than the time step used for solving.
    To avoid stiffness issues the control times will be rounded to the closest numerical solver timesteps at solving time.
    """
    function createJuMPControlVariable!(model::Model, varname::Union{Symbol,String}, cutpoints=zeros(0), values=0.0; lower_bound=-Inf, upper_bound=+Inf)
        # Preprocess varname to get both symbol and string representation
        if isa(varname,Symbol)
            varname_str = string(varname)
            varname_sym = varname
        else  #is a string
            varname_str = varname
            varname_sym = Symbol(varname)
        end

        # Check if the provided values and cutpoints are valid
        try 
            StepFunc1DData(cutpoints,values)
        catch
        end

        # If starting values is a scalar broadcast it
        if ndims(values)==0
            values = repeat([values],length(cutpoints)+1)
        end

        # Check that the variable name is not alreay registered
#=         if haskey(model, varname_sym)
            error("Variable already exists in the model").
        end =#

        # Create the corresponding JuMP variable
        # The macro only takes a directly typed expression, so I cannot use eval(varname_sym)
        # to register the variable in one command
        # https://jump.dev/JuMP.jl/stable/manual/nlp/#User-defined-functions-with-vector-inputs
        # First create an anonymous variable
        varbind = @variable(model,
            [i=1:length(cutpoints)+1], # anonymous variable with desired dimensions
            start=values[i],
            lower_bound = lower_bound,
            upper_bound = upper_bound,
            base_name=varname_str 
            )
        # Then register the correpsonding symbol in the model
        model[varname_sym] = varbind

        # Finally create the corresponding stepFuncData object and return it
        return StepFunc1DData(cutpoints,varbind)
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
        @unpack ρ,m,K,α,β,C_data = p

        # Sensitive cells dynamics
        du[1] = s′ = ρ*s*(1-((s+m*r)/K)) - α*stepFuncVal(C_data,t)*s
        # Resistant cells dynamics
        du[2] = r′ = ρ*r*(1-((s+m*r)/K)) - β*r*s/K
    end

    struct DSSolution
        u::AbstractArray#{T} where T<:Number
        du::AbstractArray#{T2} where T2<:Number
        t::AbstractArray#{T3} where T3<:Number
        p::NamedTuple
        #dynamicsFuncName::string
    end

    function DSSolution(solution::DifferentialEquations.ODESolution)::DSSolution
        # TODO
    end 

    # TODO create plot methods, compare functions

end



