function F(x, y; λ_1= 0.0, λ_2 = 0.01)
        mean_y = mean(y.instance_values, dims=2)[:,1]
        not_solved = .!y.solved_instances

        m = mean(mean_y[ not_solved ])
        if isnan(m)
            m = 0.0
        end

        m + λ_1*sum(not_solved) + λ_2*norm(x,1)
end

function f(x, y)
    length(y.solved_instances) - Float64(sum(y.evaluated_instances .* y.solved_instances))
end
