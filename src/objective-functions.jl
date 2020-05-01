function F(x, y)
        ids = y[:ids]
        mean_y = mean(y[:y], dims=2)[:,1]
        I = .!ids

        m = mean(mean_y[ I ])
        if isnan(m)
            m = 0.0
        end

        10m + 1e2sum(I) + 0.1norm(x,1)
end

function f(x, y)
    Float64(sum(y[:ids]))
end
