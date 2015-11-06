println("Loading Packages...")
using Distributions, DataFrames
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

#function rplot(xx,yy,dat,out,plotter)
#  writedlm(dat, [xx yy])
#  run(`Rscript $plotter $dat $out`)
#end

#= dp Example:
  include("dp.jl")

  @time dp(10, 1, x -> cdf(Gamma(1,1), x) , [.1,3], 11)
=#
