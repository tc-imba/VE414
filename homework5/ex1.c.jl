#=
homework5:
- Julia version: 1.1.1
- Author: liu
- Date: 2019-07-02
=#
using Distributions

function importance_sampling_approximation(x, n)
    y = rand(Normal(x, 1), n)
    w = sqrt(2 * pi) ./ (1 .+ y.^2)
    return mean(y.*w) / mean(w)
end

grid_sizes = [50, 250, 750, 1500, 3000]
for n in grid_sizes
    println(importance_sampling_approximation(0.5, n))
end
