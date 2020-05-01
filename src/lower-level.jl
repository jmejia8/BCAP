function lower_level_optimizer(Φ,problem,status,information,options,t; parameters=nothing)

    if isnothing(parameters)
        @error "Please provide a targetAlgorithm and some instances."
        return
    end



    for i = 1:length(Φ)
        if parameters.parms_type[i] <: Integer
            Φ[i] = round(Integer, Φ[i])
        elseif parameters.parms_type[i] <: AbstractFloat
            Φ[i] = round(Φ[i], digits=parameters.significant_digits)
        end
    end

    # when initialization
    if t == 0
        # seed = abs(rand(Int))
        y = call_target_algorithm(parameters.targetAlgorithm, Φ, parameters.benchmark, seed = parameters.seed, calls_per_instance=parameters.calls_per_instance)
        mean_y = mean( y, dims=2 )[:,1]
        ids = map( w -> w ≈ 0.0 ,  mean_y)

        ll_y = Dict( :y => y, :ids => ids, :ids_eval => ones(Bool,length(ids)), :feasible => true, :seed => parameters.seed )

        return Bilevel.LLResult(ll_y, problem.f(Φ, ll_y); f_calls=length(ids)*parameters.calls_per_instance)

    end

    # distances
    dists = map( sol -> norm(Φ - sol.x, 1), status.population )

    I = sortperm(dists)

    if dists[I[1]] ≈ 0.0
        options.debug && @info("Found Φ_new already in P")
        ll_y = status.population[I[1]].y

        return Bilevel.LLResult(deepcopy(ll_y), status.population[I[1]].f; f_calls = 0.0)

    end

    K = parameters.K

    y = zeros(size(status.population[I[1]].y[:y]))
    m = exp.(- dists[I[1:K]] )

    for i = 1:K
        ids_valuated = status.population[I[i]].y[:ids_eval]
        y[ids_valuated, :] += m[i]*status.population[I[i]].y[:y][ids_valuated,:]
    end

    y = y ./ sum(m)
    y = map( a -> a ≈ 0.0 ? 0.0 : a , y )

    mean_y = mean( y, dims=2 )[:,1]
    ids_to_eval = map( w -> w > 0.0,  mean_y)

    if length(ids_to_eval) > 1 && sum(ids_to_eval) > length(ids_to_eval) / 2
        I = sortperm(mean_y)
        ids_to_eval[I[1: (length(I) ÷ 2) ]] .= 0
        options.debug && print("LL approx: ")
    end


    seed = abs(rand(Int))

    options.debug && println("Solving ", sum(ids_to_eval), " instances.")
    y_new = call_target_algorithm(parameters.targetAlgorithm, Φ, parameters.benchmark, ids=ids_to_eval, seed = seed, calls_per_instance=parameters.calls_per_instance)
    f_calls = sum(ids_to_eval) * parameters.calls_per_instance
    y[ids_to_eval,:] = y_new[ids_to_eval,:]

    feasible =  sum(ids_to_eval) == length(ids_to_eval)

    mean_y = mean( y, dims=2 )[:,1]
    ids = map( w -> w ≈ 0.0 ,  mean_y)
    ll_y = Dict( :y => y, :ids => ids, :ids_eval => ids_to_eval, :feasible => feasible, :seed => seed )



    Bilevel.LLResult(ll_y, problem.f(Φ, ll_y) ; f_calls=f_calls)

end

function lower_level_optimizer(sol::Bilevel.xf_indiv, problem,status,information,options,t; parameters = nothing)
    if isnothing(parameters)
        @error "Please provide a targetAlgorithm and some instances."
        return
    end

    Φ = sol.x

    ids_to_eval = .!sol.y[:ids_eval]

    options.debug && println("Reevaluation: Solving ", sum(ids_to_eval), " instances.")


    if sum(ids_to_eval) == 0
        sol.y[:feasible] = true
        return Bilevel.LLResult(sol.y, sol.f ; f_calls=0)
    end

    y = copy(sol.y[:y])


    seed = sol.y[:seed]
    y_new = call_target_algorithm(parameters.targetAlgorithm, Φ, parameters.benchmark, ids=ids_to_eval, seed = seed, calls_per_instance=parameters.calls_per_instance)
    f_calls = sum(ids_to_eval)*parameters.calls_per_instance


    y[ids_to_eval,:] = y_new[ids_to_eval,:]


    mean_y  = mean( y, dims=2 )[:,1]
    ids_new = map( w -> w ≈ 0.0 ,  mean_y)


    ll_y = Dict( :y => y, :ids => ids_new, :ids_eval => ones(Bool, length(ids_to_eval)), :feasible => true, :seed => seed )



    Bilevel.LLResult(ll_y, problem.f(Φ, ll_y) ; f_calls=f_calls)

end


function call_target_algorithm(targetAlgorithm, Φ, benchmark; ids=ones(Bool, length(benchmark)), seed = 0, calls_per_instance = 1)

    Errors_shared = SharedArray{Float64}(length(benchmark), calls_per_instance)

    the_instances = findall(ids)
    old_seed = abs(rand(Int))
    seed!(seed)

    sd = abs.(rand(Int, calls_per_instance))
    if calls_per_instance == 1
        sd[1] = seed
    end
    @sync @distributed for r = 1:calls_per_instance
        for i = the_instances
            err = targetAlgorithm(Φ, benchmark[i], sd[r])

            Errors_shared[i, r] = err
        end
    end

    seed!(old_seed)

    Errors = Matrix(Errors_shared)

    return  Errors

end
