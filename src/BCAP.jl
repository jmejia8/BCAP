module BCAP

using Bilevel
import Statistics: mean, var, median
import LinearAlgebra: norm
import Printf: @printf
import Random: randperm
using Distributed, SharedArrays
using BiApprox, MLKernels, Optim


@everywhere import Random: seed!


include("structures.jl")
include("distributed.jl")
include("objective-functions.jl")
include("initialize.jl")
include("lower-level.jl")
include("upper-level.jl")
include("display.jl")

"""
    configure(target_algorithm::Function,
                        parameters_info::Parameters,
                        benchmark::Benchmark;
                        ll_func = f,
                        ul_func = F,
                        bcap_config = BCAP_config(),
                        debug = false,
                        store_convergence=false,
                        budget=500)

Configure a `target_algorithm` using BCAP
"""

function configure(target_algorithm::Function,
                    parameters_info::Parameters,
                    benchmark::Benchmark;
                    ll_func = f,
                    ul_func = F,
                    bcap_config = BCAP_config(),
                    debug = false,
                    store_convergence=false,
                    budget=2000)

    bounds, parameters_types = parameters_info.bounds, parameters_info.types

    bb = bounds[:, map( t -> t <: Real, parameters_types )]

    s = prod(bb[2,:] - bb[1,:] .+ 1)
    s = round(Int, 0.1s)


    D_ = size(bounds, 2)

    K = bcap_config.K
    if bcap_config.N < 1
        bcap_config.N = max(2K, min(s , K*D_))
        debug && @info("Using default population size: $(bcap_config.N)")
    end
    bcap_config.K_ll = min(bcap_config.N, max(1, bcap_config.K_ll) )
    bcap_config.parms_type = parameters_types
    bcap_config.benchmark = benchmark
    bcap_config.η_max = 1.2
    bcap_config.targetAlgorithm = target_algorithm


    options = Bilevel.Options(F_calls_limit=Inf,
                        f_calls_limit=budget,
                        F_tol=1e-5,
                        f_tol=1e-5,
                        store_convergence=store_convergence,
                        debug=debug)

    information = Bilevel.Information(f_optimum=0.0, F_optimum=0.0)

    Errors_shared = SharedArray{Float64}(length(benchmark), bcap_config.calls_per_instance)
    LL_optimizer(Φ,problem,status,information,options,t) = lower_level_optimizer(Φ,problem,status,information,options,t; parameters = bcap_config, Errors_shared = Errors_shared)




    bcap_config.p = nprocs()
    debug && @info "working with $(nprocs()) processes"


    debug && @info("Running BCAP...")
    method = Algorithm(bcap_config;
                initialize! = initialize!,
                update_state! = update_state!,
                final_stage! = final_stage!,
                lower_level_optimizer = LL_optimizer,
                is_better = is_better,
                stop_criteria = stop_criteria,
                information = information,
                options = options)




    results = Bilevel.optimize(ul_func, ll_func, bounds, bounds, method)
end


export configure, Parameters, Instance, Benchmark, BCAP_config

end # module
