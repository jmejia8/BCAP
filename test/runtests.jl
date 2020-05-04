using BCAP
using Test

function test1()
    # target_algorithm, parameters_info, instances

    target_algorithm(Φ, instance, seed = 0) = begin
        #return sum( (Φ .- (1:length(Φ))).^2 )

        if instance.index % 3 == 0
            v = rand()
        else
            v = (instance.index-1)*(sum( abs.(Φ .- (0:length(Φ)-1) )  ) + 0.5rand())
        end

        return v <= instance.optimum ? 0.0 : v
    end

    bounds = Array([ zeros(5) 10ones(5) ]')
    bounds[:, 1] = [0, 1.0]
    bounds[:, 3] = [0, 4.3]
    bounds[:, end-2] = [1, 500.0]
    parameters = Parameters(bounds, [Bool, Int, Float64, Int, Float64 ])

    benchmark = [Instance(0.5i, nothing, i) for i = 1:10]

    res = configure(target_algorithm, parameters, benchmark, debug = true )
    display(res)
    println("============================")
    for sol in res.population
        @show sol.x
        @show sol.F
        println("-----")
    end
    res.best_sol.F >= 0.0
end # function


function test2()
    # target_algorithm, parameters_info, instances

    target_algorithm(Φ, instance, seed = 0) = begin
        #return sum( (Φ .- (1:length(Φ))).^2 )

        if instance.index % 5 == 0
            v = rand()
        else
            v = (instance.index-1)*(sum( abs.(Φ .- ((1:length(Φ)) .% 2) )  ) + 0.5rand())
        end

        return v <= instance.optimum ? 0.0 : v
    end

    D = 10
    bounds = Array([ zeros(D) ones(D) ]')
    parameters = Parameters(bounds, repeat([Bool], D))

    benchmark = [Instance(0.5i, nothing, i) for i = 1:10]

    res = configure(target_algorithm, parameters, benchmark, debug = true )
    display(res)
    println("============================")
    for sol in res.population
        @show sol.x
        @show sol.F
        println("-----")
    end
    res.best_sol.F >= 0.0
end # function

@test test2()
