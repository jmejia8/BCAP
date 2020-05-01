module BCAP

using Bilevel
import Statistics: mean, var, median
import LinearAlgebra: norm
import Printf: @printf
import Random: randperm
using Distributed, SharedArrays
@everywhere import Random: seed!


include("structures.jl")
include("distributed.jl")
include("objective-functions.jl")
include("initialize.jl")
include("lower-level.jl")
include("upper-level.jl")

"""
    configure(target_algorithm, parameters_info, instances)

Configure a `target_algorithm` using BCAP
"""
function configure(target_algorithm::Function,
                    parameters_info::Parameters,
                    benchmark::Benchmark;
                    ll_func = f,
                    ul_func = F,
                    bcap_config = BCAP_config(),
                    debug = false,
                    budget=200)

    bounds, parameters_types = parameters_info.bounds, parameters_info.types
    D_ = size(bounds, 2)
    K = 6

    bcap_config.N = K*D_
    bcap_config.K = K
    bcap_config.parms_type = parameters_types
    bcap_config.significant_digits = 6
    bcap_config.calls_per_instance = 1
    bcap_config.benchmark = benchmark
    bcap_config.seed = 1
    bcap_config.η_max = 1.2
    bcap_config.targetAlgorithm = target_algorithm
    bcap_config.p = 1

    options = Bilevel.Options(F_calls_limit=Inf,
                        f_calls_limit=budget*length(benchmark),
                        F_tol=1e-5,
                        f_tol=1e-5,
                        store_convergence=false,
                        debug=debug)

    information = Bilevel.Information(f_optimum=0.0, F_optimum=0.0)

    LL_optimizer(Φ,problem,status,information,options,t) = lower_level_optimizer(Φ,problem,status,information,options,t; parameters = bcap_config)

    bcap_config.p > 1 && addprocs(bcap_config.p - 1)
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


export configure, Parameters, Instance, Benchmark

end # module
