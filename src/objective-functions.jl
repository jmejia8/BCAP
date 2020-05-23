function F(x, y; λ=0.01)
    mean_y = mean(y.instance_values, dims=2)[:,1]
    not_solved = .!y.solved_instances
    r = λ*norm(x,1)

    if sum(not_solved) == 0
        return r
    end

    m = (mean_y[ not_solved ])

    I = m .> 1000

    m[I] = min.(900, 10log10.(m[I]))


    mean(m) + 1e3sum(not_solved) + r
end

function f(x, y)
    length(y.solved_instances) - Float64(sum(y.evaluated_instances .* y.solved_instances))
end
