function lower_level_optimizer(
    Φ,
    problem,
    status,
    information,
    options,
    t;
    parameters = nothing,
    Errors_shared = SharedArray{Float64}(length(benchmark), 31),
)

    if isnothing(parameters)
        @error "Please provide a targetAlgorithm and some instances."
        return
    end



    # convert types of parameters
    for i = 1:length(Φ)
        if parameters.parms_type[i] <: Integer
            Φ[i] = round(Integer, Φ[i])
        elseif parameters.significant_digits > 0 &&
               parameters.parms_type[i] <: AbstractFloat
            Φ[i] = round(Φ[i], digits = parameters.significant_digits)
        end
    end

    num_instances = length(parameters.benchmark)

    # initialization step
    if t == 0
        # seed = abs(rand(Int))
        y = call_target_algorithm(
            parameters.targetAlgorithm,
            Φ,
            parameters.benchmark,
            seed = parameters.seed,
            calls_per_instance = parameters.calls_per_instance,
            Errors_shared = Errors_shared,
        )

        ll_y = LLSolution(
            instance_values = y,
            evaluated_instances = ones(Bool, num_instances),
            seed = parameters.seed,
        )

        return Bilevel.LLResult(
            ll_y,
            problem.f(Φ, ll_y);
            f_calls = num_instances * parameters.calls_per_instance,
        )

    end

    # distances { ‖ Φ_new - Φ ‖ Φ ∈ P }
    distances = map(sol -> norm(Φ - sol.x, 1), status.population)

    I = sortperm(distances)

    if distances[I[1]] ≈ 0.0
        options.debug && @info("Found Φ_new already in P")
        ll_y = status.population[I[1]].y

        return Bilevel.LLResult(
            deepcopy(ll_y),
            status.population[I[1]].f;
            f_calls = 0.0,
        )

    end

    K = parameters.K_ll

    y = zeros(size(status.population[I[1]].y.instance_values))
    m = exp.(-distances[I[1:K]])

    for i = 1:K
        valuated_instances = status.population[I[i]].y.evaluated_instances
        y[valuated_instances, :] +=
            m[i] *
            status.population[I[i]].y.instance_values[valuated_instances, :]
    end

    y = y ./ sum(m)
    y = map(a -> a ≈ 0.0 ? 0.0 : a, y)

    # display(y)

    mean_y = mean(y, dims = 2)[:, 1]
    instances2eval = map(w -> w > 0.0, mean_y)

    if length(instances2eval) > 1 &&
       sum(instances2eval) > length(instances2eval) / 2
        I = sortperm(mean_y)
        instances2eval[I[1:(length(I)÷2)]] .= 0
        options.debug && println("Φ = ", Φ, " is too bad")
    end


    options.debug && println("Solving ", sum(instances2eval), " instances.")

    seed = parameters.seed

    y_new = call_target_algorithm(
        parameters.targetAlgorithm,
        Φ,
        parameters.benchmark,
        ids = instances2eval,
        seed = seed,
        calls_per_instance = parameters.calls_per_instance,
        Errors_shared = Errors_shared,

    )


    f_calls = sum(instances2eval) * parameters.calls_per_instance
    y[instances2eval, :] = y_new[instances2eval, :]

    ll_y = LLSolution(
        instance_values = y,
        evaluated_instances = instances2eval,
        seed = seed,
    )

    Bilevel.LLResult(ll_y, problem.f(Φ, ll_y); f_calls = f_calls)

end

function lower_level_optimizer(
    sol::Bilevel.xf_indiv,
    problem,
    status,
    information,
    options,
    t;
    parameters = nothing,
    Errors_shared = nothing,
)
    if isnothing(parameters)
        @error "Please provide a targetAlgorithm and some instances."
        return
    end

    Φ = sol.x
    lower_level = sol.y
    instances2eval = .!lower_level.evaluated_instances

    options.debug &&
        println("Reevaluation: Solving ", sum(instances2eval), " instances.")


    if sum(instances2eval) == 0
        lower_level.isfeasible = true
        return Bilevel.LLResult(lower_level, sol.f; f_calls = 0)
    end

    y = copy(lower_level.instance_values)

    seed = lower_level.seed
    y_new = call_target_algorithm(
        parameters.targetAlgorithm,
        Φ,
        parameters.benchmark,
        ids = instances2eval,
        seed = seed,
        calls_per_instance = parameters.calls_per_instance,
        Errors_shared = Errors_shared,
    )
    f_calls = sum(instances2eval) * parameters.calls_per_instance


    y[instances2eval, :] = y_new[instances2eval, :]

    ll_y = LLSolution(
        instance_values = y,
        evaluated_instances = ones(Bool, length(instances2eval)),
        seed = seed,
    )

    Bilevel.LLResult(ll_y, problem.f(Φ, ll_y); f_calls = f_calls)

end



function call_target_algorithm_parallel(
    targetAlgorithm,
    Φ,
    benchmark;
    ids = ones(Bool, length(benchmark)),
    seed = 1,
    calls_per_instance = 1,
    Errors_shared = nothing,
)





    the_instances = findall(ids)
    old_seed = abs(rand(Int))

    sd = seed * ones(Int, calls_per_instance)
    sd += 0:(length(sd)-1)

    L = [(r, i) for r = 1:calls_per_instance for i in the_instances]

    @sync @distributed for a in L
        r, i = a
        err = targetAlgorithm(Φ, benchmark[i], sd[r])

        Errors_shared[i, r] = err

    end

    seed!(old_seed)

    Errors = zeros(length(benchmark), calls_per_instance)
    Errors[the_instances, :] = Errors_shared.s[the_instances, :]


    return Errors

end

function call_target_algorithm(
    targetAlgorithm,
    Φ,
    benchmark;
    ids = ones(Bool, length(benchmark)),
    seed = 1,
    calls_per_instance = 1,
    Errors_shared = nothing,
)


    if nprocs() > 1
        return call_target_algorithm_parallel(
            targetAlgorithm,
            Φ,
            benchmark;
            ids = ids,
            seed = seed,
            calls_per_instance = calls_per_instance,
            Errors_shared = Errors_shared,
        )
    end

    the_instances = findall(ids)
    Errors = zeros(length(benchmark), calls_per_instance)

    old_seed = abs(rand(Int))
    seed!(seed)



    sd = seed * ones(Int, calls_per_instance)
    sd += 0:(length(sd)-1)
    for r = 1:calls_per_instance
        for i in the_instances
            err = targetAlgorithm(Φ, benchmark[i], sd[r])

            Errors[i, r] = err
        end
    end

    seed!(old_seed)


    return Errors

end
