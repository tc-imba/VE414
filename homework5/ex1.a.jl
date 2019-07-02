#=
homework5:
- Julia version: 1.1.1
- Author: liu
- Date: 2019-07-02
=#
using Distributions

function reject_sampling_approximation(x, n)
    A = Array{Float64}(undef, n)
    i = 0
    while i < n
        v = rand(1)[1]
        y = rand(Normal(x, 1), 1)[1]
        if v <= 1 / (1 + y^2)
            A[i+=1] = y
        end
    end
    return mean(A)
end

grid_sizes = [50, 250, 750, 1500, 3000]
for n in grid_sizes
    println(reject_sampling_approximation(0.5, n))
end
