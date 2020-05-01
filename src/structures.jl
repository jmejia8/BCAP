mutable struct Parameters
    bounds::Matrix{Float64}
    types::Vector
end # mutable struct

mutable struct Instance
    optimum::Float64
    value
    index::Int

end # mutable struct

Benchmark = Array{Instance}


mutable struct BCAP_config
    N::Int64
    K::Int64
    η_max::Float64
    λ::Float64
    parms_type::Array
    significant_digits::Int
    calls_per_instance::Int
    targetAlgorithm::Function
    benchmark::Array{Instance}
    seed::Int
    training_population::Set
    approx_model
    F_approx::Function
    max_global_iters::Int
    max_local_iters::Int
    p::Int
end

function BCAP_config(; N = 30,
                K = 6,
                η_max = 1.0,
                λ = 1e-1,
                parms_type = [],
                significant_digits = 6,
                calls_per_instance = 1,
                targetAlgorithm = identity,
                benchmark = Instance[],
                seed = 0,
                training_population = Set([]),
                approx_model = nothing,
                F_approx = identity,
                max_global_iters = 100,
                p = 1,
                max_local_iters = 10)

    BCAP_config(N,
                K,
                η_max,
                λ,
                parms_type,
                significant_digits,
                calls_per_instance,
                targetAlgorithm,
                benchmark,
                seed,
                training_population,
                approx_model,
                F_approx,
                max_global_iters,
                max_local_iters,
                p)
end
