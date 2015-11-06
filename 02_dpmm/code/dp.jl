println("Loading Packages...")
using Gadfly, Distributions, DataFrames
println("Completed!")

run(`mkdir -p temp`)

function rdir(N,a)
  K = size(a)
  if K == ()
    K = 1
  else
    K = K[1]
  end
  x = zeros(N,K)
  for i in 1:N, j in 1:K
    x[i,j] = rand(Gamma(a[j],1)) # shape and scale
  end
  rowsums = x * ones(K)
  out = x ./ rowsums

  return out
end

function dp(N,a,pG,xlim,J=100) # if discrete, set J = size(xlim)
  x = linspace(xlim[1],xlim[2],J)
  dG0 = [ pG(x[1]); Any[pG(x[j]) - pG(x[j-1]) for j in 2:J] ]
  out = rdir(N, a .* dG0)
  G = cumsum(out')' # cumsum only works my column.
  return (G,x)
end

function rplot(x,y)
  file = "temp/temp.dat"
  writedlm(file, [x y])
  run(`Rscript plot.R $file`)
end

#= dp Example:
  include("dp.jl")
  pG(x) = cdf(Gamma(1,1),x)
  xlim = [.1,3]
  N,a = 3,1
  dp(N,a,pG,xlim)

  # Plotting:
  plot(x=[1,2,3], y=[2,3,4])
  draw(PDF("./test.pdf"), Gadfly.plot(x=[1,2,3],y=[1,2,3]))
=#
