function F(x, y)
        mean_y = mean(y.instance_values, dims=2)[:,1]
        not_solved = .!y.solved_instances

        m = norm(mean_y[ not_solved ])
        if isnan(m)
            m = 0.0
        end

        m + 0.1norm(x,1)
end

function f(x, y)
    Float64(sum(y.solved_instances))
end
