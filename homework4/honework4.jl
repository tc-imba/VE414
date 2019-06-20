function grid_approximation(x, n)
  y_grid = collect(range(-pi/2 , length=n , stop=pi/2))
  unnormalised_posterior = map(u->exp(-(x-tan(u))^2/2), y_grid)
  unnormalised_expectation = map(u->exp(-(x-tan(u))^2/2)*tan(u), y_grid)
  A = pi * sum(unnormalised_posterior) / n
  E = pi * sum(unnormalised_expectation) / A / n
  return A, E
end
grid_sizes = [50, 250, 750, 1500, 3000]
for n in grid_sizes
  println(grid_approximation(0.5, n))
end

using Distributions
function direct_grid_approximation(x, n, m)
  y_grid = collect(range(-pi/2 , length=n , stop=pi/2))
  unnormalised_posterior = map(y->exp(-(x-tan(y))^2/2), y_grid)
  A = pi * sum(unnormalised_posterior) / n
  posterior = unnormalised_posterior / A
  samples = map(y->tan(y), wsample(y_grid, posterior, m))
  E = sum(samples) / m
  return A, E
end
grid_sizes = [50, 250, 750, 1500, 3000]
sample_sizes = [100, 1000]
for n in grid_sizes
  for m in sample_sizes
    println(direct_grid_approximation(0.5, n, m))
  end
end

function no_transform_grid_approximation(x, a, b, n)
  if n <= 1000
    y_grid = collect(range(a, length=n, stop=b))
    newa = a
  elseif n > 1000 && n <= 2000
    nm = 1000
    na = round(Int, (n-nm)/2)
    l = (b-a)/(nm-1)
    newa = a - l*na
    y_grid = collect(range(newa, step=l, length=n))
  else n > 2000
    nm = round(Int, n/2)
    na = round(Int, (n-nm)/2)
    l = (b-a)/(nm-1)
    newa = a - l*na
    y_grid = collect(range(newa, step=l, length=n))
  end
  unnormalised_posterior = map(y->exp(-(x-y)^2/2)/(1+y^2), y_grid)
  unnormalised_expectation = map(y->exp(-(x-y)^2/2)/(1+y^2)*y, y_grid)
  A = (y_grid[n]-y_grid[1]) * sum(unnormalised_posterior) / n
  E = (y_grid[n]-y_grid[1]) * sum(unnormalised_expectation) / A / n
  return A, E
end

grid_sizes = [50, 250, 750, 1500, 3000]
for n in grid_sizes
  println(no_transform_grid_approximation(0.5, -5, 5, n))
end
