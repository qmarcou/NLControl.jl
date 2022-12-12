"""
Collection of functions implementing ODE dynamics using DifferentialEquations.jl style.
"""
module ODESystems
    export twoPopRSC!
    # check this package out: https://juliapackages.com/p/parameterizedfunctions
    using Parameters
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

end