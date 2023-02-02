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
    using JuMP: Model,@variable,@variables,fix,@NLexpressions,@NLconstraint,add_nonlinear_expression,NonlinearExpression,start_value,register
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
        # constraints
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
        function dynamics(u,p,t)
            du = similar(u)
            system.dynamics!(du,u,p,t)
            return du
        end
        return dynamics
    end

    function generateSplattedCopyDynamics(system::DynamicalSystem,
        controlVarNames::AbstractArray{<:Union{String,Symbol}})::Function
        dynamics = generateCopyDynamics(system)
        if length(controlVarNames)!=1
            error("Handling of multiple control variables is not implemented yet.")
        end
        ndims = system.ndims
        p_copy = copy(system.p)
        function splattedDynamics(t,x...)
            # assign current control value to p_copy (need to convert to array using collect, x being a tuple)
            p_copy[controlVarNames[1]] = StepFunc1DData(p_copy[controlVarNames[1]].cutpoints,collect(x[ndims+1:end]))
            return dynamics(collect(x[1:ndims]), # u
                    p_copy, # p with updated control term
                    t)
        end
        return splattedDynamics
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

    """
    memoize(func::Function, n_outputs::Int)

    Take a function `func` and return a vector of length `n_outputs`, where element
    `i` is a function that returns the equivalent of `func(x...)[i]`.

    To avoid duplication of work, cache the most-recent evaluations of `func`.
    Because `func_i` is auto-differentiated with ForwardDiff, our cache needs to
    work when `x` is a `Float64` and a `ForwardDiff.Dual`.
    """
#=     function memoizeDynamics(splattedDyn::Function, n_outputs::Int)
        last_x, last_f, last_t = nothing, nothing, nothing
        last_dx, last_dfdx, last_dt = nothing, nothing, nothing
        function splattedDyn_i(i,p,t,x::T...) where {T<:Real}
            if T == Float64
                if (x != last_x) || (t != last_t)
                    last_x, last_f, last_t = x, splattedDyn(x...), t
                end
                return last_f[i]::T
            else
                if (x != last_dx) || (t != last_dt) || (!isa(last_dfdx,T)) 
                    # the last check seems necessary to prevent Forward diff from crashing
                    last_dx, last_dfdx, last_dt = x, splattedDyn(x...), t
                end
                return last_dfdx[i]::T
            end
        end
        return [(p,t,x...) -> splattedDyn_i(i,p,t,x...) for i in 1:n_outputs]
    end =#
    function memoize(func::Function, n_outputs::Int)
        last_x, last_f = nothing, nothing
        last_dx, last_dfdx = nothing, nothing
        function func_i(i,x::T...) where {T<:Real}
            if T == Float64
                if (x != last_x)
                    last_x, last_f = x, func(x...)
                end
                return last_f[i]::T
            else
                if (x != last_dx) || (!isa(last_dfdx,T)) 
                    # the last check seems necessary to prevent Forward diff from crashing
                    last_dx, last_dfdx = x, func(x...)
                end
                return last_dfdx[i]::T
            end
        end
        return [(x...) -> func_i(i,x...) for i in 1:n_outputs]
    end


    """

    Create a base JuMP model based on the provided dynamical system.
    The returned model object only contains dynamics and starting solutions
    """
    function addJuMPNLDynamics!(model::Model,
        system::DynamicalSystem,
        Δt,
        u0,
        controlVarNames::AbstractArray{<:Union{String,Symbol}})

        # Get array sizes
        n_rows = system.ndims
        n_cols = length(Δt) + 1
        
        @assert length(u0) == n_rows

        # Make a working copy of the parameters Dict
        p_copy = copy(system.p)

        # Extract starting values from JuMP control variables
        if length(controlVarNames)>1
            error("Handling of multiple control variables is not implemented yet.")
        end
        for var in controlVarNames
            if isa(var,String)
                # DynamicalSystem.p is a Dict{Symbol,Any}
                var = Symbol(var)
            end
            symbolic_val = p_copy[var]
            p_copy[var] = StepFunc1DData(symbolic_val.cutpoints,start_value.(symbolic_val.values))
        end

        # Make a copy of the system containing the parameters copy
        system_copy = DynamicalSystem(p_copy,system.dynamics!,system.ndims)

        # Compute an initializing solution using the control starting values
        u_init = solveEuler(system_copy,u0,Δt).u

        # Create the variables and initialize them 
        @variables(model,begin
            u[i=1:n_rows,j=1:n_cols], (start=u_init[i,j])
            t[j=1:n_cols] # create a time vector variable 
        end)
        
        # Fix the vector time variable
        tmp_t = append!([0.0],cumsum(Δt))
        for j in 1:n_cols
            fix(t[j], tmp_t[j])
        end

        # Fix initial state to u0
        for i in 1:n_rows
            fix(u[i,1], u0[i]; force = true) 
        end

        # Create a Matrix of NLexpression containing the dynamics 
        # We have both vector input and vector output functionsthat we have to circumvent
        # Vector input: https://jump.dev/JuMP.jl/stable/manual/nlp/#User-defined-functions-with-vector-inputs
        # Vector output: https://jump.dev/JuMP.jl/stable/manual/nlp/#More-complicated-examples
        splattedDyn = generateSplattedCopyDynamics(system,controlVarNames)
        memoizedDyn = memoize(splattedDyn, system.ndims)
        du_binding = Matrix{NonlinearExpression}(undef,n_rows,n_cols)
        for i in 1:system.ndims 
            register(model, Symbol("du_"*string(i)),
             1+system.ndims+length(system.p[controlVarNames[1]].values), # time + system + control variables dimensions
              memoizedDyn[i]; autodiff = true)
            for j in 1:n_cols
                expr = :($(Symbol("du_"*string(i)))($(t[j]),$(u[:,j])...,$(system.p[controlVarNames[1]].values)...))
                du_binding[i,j] = add_nonlinear_expression(model,expr)
            end
        end
        # register the matrix as a variable in the model
        model[:du] = du_binding

        # Finally add the integration constraint using the expression of the dynamics
        register(model,:eulerStep,3,NumInt.eulerStep,autodiff=true)
        @NLconstraint(model, numint[i=1:n_rows,j=2:n_cols], u[i,j] == eulerStep(u[i,j - 1], du_binding[i,j - 1], t[j]-t[j-1]))

        return (model,u,du_binding,t)
    end

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
        if haskey(model, varname_sym)
            error("Variable '"*varname_str*"' already exists in the model.")
        end

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



