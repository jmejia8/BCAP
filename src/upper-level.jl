function update_state!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    t_main_loop,
)

    if status.stop
        return
    end

    best_feasible = deepcopy(status.best_sol)
    Bilevel.BCAOperators.update_state!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        t_main_loop,
    )

    force_reevaluation = sum(status.best_sol.y.solved_instances) == length(status.best_sol.y.solved_instances)
    options.debug && force_reevaluation && @info "Re-evualing since infeasible elite solution seems true optimum."

    if !status.best_sol.y.isfeasible || status.best_sol.F > best_feasible.F
        status.best_sol = best_feasible
    end

    if status.best_sol.y.isfeasible && status.best_sol.f == 0.0
        status.stop = true
        status.stop_msg = "Optimum found 1"
        return
    end

    reevaluate = status.stop || (t_main_loop > 0 && t_main_loop % parameters.t_reevaluation == 0)
    options.debug && reevaluate && @info "Re-evualing since ($(t_main_loop)th iter) % $(parameters.t_reevaluation) = 0"

    # force reevaluation of last population
    last_iteration = options.f_calls_limit - status.f_calls <= (parameters.N) *
                     length(parameters.benchmark) *
                     parameters.calls_per_instance

    force_reevaluation = force_reevaluation || last_iteration
    options.debug && force_reevaluation && @info "Re-evualing since last iteration."

    if !reevaluate && !force_reevaluation
        status.final_time = time()
        return
    end

    # replace infeasible
    insert_best_to_pop = true
    id_worst = 1
    for i in 1:length(status.population)
        sol = status.population[i]
        if !sol.y.isfeasible
            p = sol.x

            pred = sol.y.solved_instances

            ll_result = engine.lower_level_optimizer(
                sol,
                problem,
                status,
                information,
                options,
                t_main_loop,
            )
            status.f_calls += ll_result.f_calls
            q = ll_result.y
            FF = problem.F(p, q)
            status.F_calls += 1

            sol.y = q
            sol.F = FF
            sol.f = ll_result.f

            push!(parameters.solutions, deepcopy(sol))
        end

        if is_better(sol, status.best_sol)
            options.debug && @info "Best sol. updated in reevaluation"
            status.best_sol.F = sol.F
            status.best_sol.f = sol.f
            status.best_sol.x = copy(sol.x)
            status.best_sol.y = deepcopy(sol.y)
            insert_best_to_pop = false
        end

        if is_better(status.population[id_worst], sol)
            id_worst = i
        end

        if status.f_calls >= options.f_calls_limit || status.best_sol.f == 0.0
            status.stop = true
            if status.best_sol.f == 0.0
                status.stop_msg = "Optimum found 2"
            else
                status.stop_msg = "Budget spent"
            end
            break
        end
    end

    if insert_best_to_pop
        status.population[id_worst] = status.best_sol
    end


    ##################################################
    ##################################################
    ##################################################
    ##################################################

    parameters.surrogated && !force_reevaluation && surrogate!(problem,
        engine,
        parameters,
        status,
        information,
        options,
        t_main_loop,)

    status.final_time = time()

end

function surrogate!(problem,
    engine,
    parameters,
    status,
    information,
    options,
    t_main_loop,
)

    unique!(parameters.solutions)


    a = problem.bounds_ul[1,:]
    b = problem.bounds_ul[2,:]



    n = length(parameters.solutions)
    n_train = min(1000, n)

    options.debug && @info "Training with ($(n_train) different confs.) / ($n in total)"

    sols = parameters.solutions[ randperm(n)[1:n_train] ]

    X = map(sol -> sol.x', sols)
    y = map(sol -> sol.F, sols)
    X = vcat(X...)
    X = (X .- a') ./ (b - a)'
    method = KernelInterpolation(y, X, λ = parameters.λ, kernel = PolynomialKernel())
    train!(method)
    F̂ = approximate(method)


    x_initial = (status.best_sol.x - a) ./ (b - a)
    optimizer = Optim.Fminbox(Optim.BFGS())
    res = Optim.optimize(F̂, zeros(length(a)), ones(length(a)), x_initial, optimizer, Optim.Options(outer_iterations = 1))
    p = a .+ (b - a) .* res.minimizer

    ll_result = engine.lower_level_optimizer(p,problem,status,information,options,0)

    status.f_calls += ll_result.f_calls
    q = ll_result.y
    FF = problem.F(p, q)
    if FF < status.best_sol.F
        status.best_sol.F = FF
        status.best_sol.f = ll_result.f
        status.best_sol.x = p
        status.best_sol.y = q
    else
        options.debug && @warn "Fail improvement (best)!"
    end

    i_worst = argmax(map( ind -> ind.F, status.population ))

    if FF < status.population[i_worst].F
        status.population[i_worst].F = FF
        status.population[i_worst].f = ll_result.f
        status.population[i_worst].x = p
        status.population[i_worst].y = q
    else
        options.debug && @warn "Fail improvement!"
    end

    push!(parameters.solutions, deepcopy(status.population[i_worst]))
end

is_better_approx(solution_1, solution_2) = solution_1.F < solution_2.F

function is_better(solution_1, solution_2)
    return solution_1.F < solution_2.F
end

function stop_criteria(status, information, options)
    if Bilevel.stop_check(status, information, options)
        return true
    end

    if status.best_sol.y.isfeasible && status.best_sol.f == 0.0
        status.stop = true
        status.stop_msg = "Optimum found 3"
        return true
    end

    Fs = map(sol -> sol.F, status.population)
    status.stop = var(Fs) < 1e-8
    status.stop_msg = status.stop ? "var(F(Φ_i)) ≈ 0" : ""

    status.stop


end

function final_stage!(status, information, options)
    i = 1
    ids = Int[]
    for sol in status.population
        if sol.y.isfeasible && is_better(sol, status.best_sol)
            status.best_sol = sol
        end

        if !sol.y.isfeasible
            push!(ids, i)
        end

        i += 1
    end

    options.debug &&
        length(ids) > 0 &&
        @info "Removing $(length(ids)) infeasible solution(s)."

    deleteat!(status.population, ids)
end
