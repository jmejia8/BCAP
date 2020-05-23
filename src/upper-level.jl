function update_state!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    t_main_loop,
)

    best = deepcopy(status.best_sol)
    Bilevel.BCAOperators.update_state!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        t_main_loop,
    )

    if !status.best_sol.y.isfeasible || status.best_sol.F > best.F
        status.best_sol = best
    end

    if status.best_sol.y.isfeasible && status.best_sol.f == 0.0
        status.stop = true
        status.stop_msg = "isfeasible or optimum found"
        status.stop_msg = "Optimum found"
        return
    end

    reevaluate = status.stop || (t_main_loop > 0 && t_main_loop % 5 != 0)

    # force reevaluation of last population
    force_reevaluation =
        options.f_calls_limit - status.f_calls <=
        (parameters.N) *
        length(parameters.benchmark) *
        parameters.calls_per_instance


    if reevaluate && !force_reevaluation
        status.final_time = time()
        return
    end

    options.debug && @info "Re-evualing...."
    options.debug && force_reevaluation && @info "Forced Re-evualuation...."



    for sol in status.population
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
        end

        if is_better(sol, status.best_sol)
            options.debug && @info "Best sol. updated in reevaluation"
            status.best_sol.F = sol.F
            status.best_sol.f = sol.f
            status.best_sol.x = copy(sol.x)
            status.best_sol.y = deepcopy(sol.y)
        end

        if status.f_calls >= options.f_calls_limit || status.best_sol.f == 0.0
            status.stop = true
            status.stop_msg = "f_calls limited or optimum found"
            break
        end
    end


    ##################################################
    ##################################################
    ##################################################
    ##################################################

    parameters.surrogated && surrogate!(problem,
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

    a = problem.bounds_ul[1,:]
    b = problem.bounds_ul[2,:]

    X = map(sol -> sol.x', status.population)
    y = map(sol -> sol.F, status.population)
    X = vcat(X...)
    X = (X .- a') ./ (b - a)'
    method = KernelInterpolation(y, X, λ = 1e-5, kernel = PolynomialKernel())
    train!(method)
    F̂ = approximate(method)


    x_initial = (status.best_sol.x - a) ./ (b - a)
    optimizer = Optim.Fminbox(Optim.BFGS())
    res = Optim.optimize(F̂, zeros(length(a)), ones(length(a)), x_initial, optimizer, Optim.Options(outer_iterations = 1))
    p = a .+ (b - a) .* res.minimizer

    ll_result = engine.lower_level_optimizer(p,problem,status,information,options,0)
    @show p

    status.f_calls += ll_result.f_calls
    q = ll_result.y
    FF = problem.F(p, q)
    if FF < status.best_sol.F
        status.best_sol.F = FF
        status.best_sol.f = ll_result.f
        status.best_sol.x = p
        status.best_sol.y = q
    else
        @info "Fail improvement (best)!"
    end

    i_worst = argmax(y)

    if FF < status.population[i_worst].F
        status.population[i_worst].F = FF
        status.population[i_worst].f = ll_result.f
        status.population[i_worst].x = p
        status.population[i_worst].y = q
    else
        @info "Fail improvement!"
    end
end

is_better_approx(solution_1, solution_2) = solution_1.F < solution_2.F

function is_better(solution_1, solution_2)
    return solution_1.F < solution_2.F #|| solution_1.f < solution_2.f

end

function stop_criteria(status, information, options)
    if Bilevel.stop_check(status, information, options)
        return true
    end

    if status.best_sol.y.isfeasible && status.best_sol.f == 0.0
        status.stop = true
        status.stop_msg = "Optimum found"
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
