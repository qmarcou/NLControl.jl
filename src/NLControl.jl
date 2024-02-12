module NLControl
    include("NumInt.jl")
    include("ODESystems.jl")
    import .ODESystems,.NumInt
end
