function update_state!(problem,engine,parameters,status,information,options,t_main_loop)

    best = deepcopy(status.best_sol)
    P_old = copy(status.population)
    Bilevel.BCAOperators.update_state!(problem,engine,parameters,status,information,options,t_main_loop)

    if !status.best_sol.y.isfeasible || status.best_sol.F > best.F
        status.best_sol = best
    end

    if status.best_sol.y.isfeasible && status.best_sol.f == length(parameters.benchmark)
        status.stop = true
        status.stop_msg = "isfeasible or optimum found"
        status.stop_msg = "Optimum found"
        return
    end

    reevaluate = status.stop || (t_main_loop > 0 && t_main_loop % 5 != 0)

    # force reevaluation of last population
    force_reevaluation = options.f_calls_limit - status.f_calls <= (parameters.N)*length(parameters.benchmark) * parameters.calls_per_instance


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

            ll_result = engine.lower_level_optimizer(sol,problem,status,information,options, t_main_loop)
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

        push!(P_old, sol)

        if status.f_calls >= options.f_calls_limit || status.best_sol.f == length(parameters.benchmark)
            status.stop = true
            status.stop_msg = "f_calls limited or optimum found"
            break
        end
    end


    if length(P_old) > parameters.N
        sort!(P_old, lt = (a, b) -> is_better(a, b))
        status.population = P_old[1:parameters.N]
    end

    status.final_time = time()

end


is_better_approx(solution_1, solution_2) = solution_1.F < solution_2.F

function is_better(solution_1, solution_2)
    return solution_1.F < solution_2.F || solution_1.f > solution_2.f

end

function stop_criteria(status, information, options)
    if Bilevel.stop_check(status, information, options)
        return true
    end

    Fs = map(sol -> sol.F, status.population )
    status.stop = var(Fs) < 1e-8
    status.stop_msg = status.stop ? "var(F(Φ_i)) ≈ 0" : ""

    status.stop


end

function final_stage!(status, information, options)
    i = 1
    ids = Int[]
    for sol = status.population
        if sol.y.isfeasible && is_better(sol, status.best_sol)
            status.best_sol = sol
        end

        if !sol.y.isfeasible
            push!(ids, i)
        end

        i += 1
    end

    options.debug && length(ids) > 0 && @info "Removing $(length(ids)) infeasible solution(s)."

    deleteat!(status.population, ids)
end
