"""
Collection of functions to perform numerical integration/ODE solving.

Though there probably exists much better implementations the point of this
module is to provide functions that are autodifferentiable and can be used
to build constraints on dynamics for JuMP models.
"""
module NumInt
    export eulerSolve
    import JuMP: AbstractJuMPScalar, GenericQuadExpr

    #function eulerStep!(u_next::Array{float64},u::Array{float64},du::Array{float64},Δx::Array{float64})
    function eulerStep!(u_next::Union{T,AbstractArray{T}},
                        u::Union{T,AbstractArray{T}},
                        du::Union{T,AbstractArray{T}},
                        Δx::Number) where T<: Union{Number,AbstractJuMPScalar}
        if ndims(u)==0
            u_next = u+Δx*du
        else
            u_next .= u.+Δx.*du
        end
    end

    function eulerStep(u::Union{T,AbstractArray{T}},
                        du::Union{T,AbstractArray{T}},
                        Δx::Number) where T<: Union{Number,AbstractJuMPScalar}
        if ndims(u)==0
            return u+Δx*du
        else
            return u.+Δx.*du
        end
    end

    function eulerSolve(systemFunc!::Function,
                        p,
                        u_0::AbstractArray{T},
                        Δt::AbstractArray{T}) where T<: Number
        # Assert time steps and initial condition are vector like
        @assert ndims(Δt) == 1
        @assert ndims(u_0) == 1
        # Initialize values
        u = similar(u_0,Union{Number,AbstractJuMPScalar},length(u_0),length(Δt)+1)
        u[:,1] = u_0
        du = similar(u)
        t = similar(Δt,length(Δt)+1)
        t[1] = 0.0
        i::Int = 1
        for Δ in Δt
            systemFunc!(view(du,:,i),view(u,:,i),p,t[i])
            eulerStep!(view(u,:,i+1),view(u,:,i),view(du,:,i),Δ)
            t[i+1]=t[i]+Δ
            i+=1
        end
        return (u=u,du=du,t=t)
    end

    function RK4Step(u_next,u,du,Δt)
        Nothing
    end

    function ()
        
    end

    function rect_integral(f,Δx) 
        integral = 0.0
        for (idx, val) in enumerate(f[1:lastindex(f)-1])
            integral += Δx*val
        end
        return integral
    end
    
    function trap_integral(f,Δx)
        integral = 0.0
        for (idx, val) in enumerate(f[1:lastindex(f)-1])
            integral += Δx*(val+f[idx+1])/2.0
        end
        return integral
    end
    
end