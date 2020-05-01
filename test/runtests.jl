using BCAP
using Test

function test1()
    # target_algorithm, parameters_info, instances

    target_algorithm(Φ, instance, seed = 0) = begin
        return sum(Φ) + 0.1rand()
    end

    bounds = Array([ zeros(5) 10ones(5) ]')
    parameters = Parameters(bounds, [ Float64, Int, Int, Float64, Int ])

    benchmark = [Instance(0.0, nothing, i) for i = 1:10]

    res = configure(target_algorithm, parameters, benchmark, debug = true )
    display(res)
    res.best_sol.F >= 0.0
end # function

@test test1()
