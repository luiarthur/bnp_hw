using Distributions, DataFrames

y = [rand(Normal(3,.1),30);
     rand(Normal(5,.1),16); 
     rand(Normal(6,.1),26)]
#=
  F(theta) = N(theta,.1²)
  G₀ = N(0,1)
=#

n = size(y)[1]
B = 10000
c = ones(Int64,B,n)
t = ones(B,n)
a = 1

for b in 2:B
  c[b,:] = c[b-1,:]
  t[b,:] = t[b-1,:]
  for i in 1:n
    u = unique(c)
    k = size(u)[1]
    
    #aux = Any[Dict(u[j] => ) for j in 1:n]
    
  end # for i in 1:n
end # for b in 1:B
