println("Loading Packages...")
  using Distributions, DataFrames#, UnicodePlots
  #https://github.com/Evizero/UnicodePlots.jl
println("Completed!")

function sort_table(dict::Dict{Any,Any}, prop=false)
  ord = sortperm( collect(values(dict)), rev=true )
  n = length(dict)
  sumv = sum(values(dict))
  out = [ ( collect( keys(dict) )[o], collect( values(dict))[o] / ifelse(prop,sumv,1) ) for o in ord]

  return out
end

function table(x)
  dict = Dict()

  for u in x
    if haskey(dict, u)
      dict[u] += 1
    else
      dict[u] = 1
    end
  end

  return dict
end

#y = [rand(Normal(3 ,s),30);
#     rand(Normal(1 ,s),3); 
#     rand(Normal(20,s),15)]

y = squeeze(readdlm("../../dat/hw2.dat")',1)
#histogram(y)
#=
  F(theta) = N(theta,.1²)
  G₀ = N(0,1)
=#

n = length(y)
B = 100
burn = Int64(B * .3)
t = ones(B,n)
a = 1
p = 0
s = 1

G0_kernel = Normal(0,3)
G0 = rand(G0_kernel, B)
lg0(x) = logpdf(G0_kernel, x )

f(x, theta) =     pdf( Normal(theta,s), x )
lf(x, theta) = logpdf( Normal(theta,s), x )

acc_t = 0
cs = 10

for b in 2:B
  t[b,:] = t[b-1,:]

  # Update θ_i | θ_{-i}
  for i in 1:n
    tb = t[b,:]
    ti = tb[i]
    #ut_wo_ti = unique( tb[find(x-> x != i, 1:n)] )
    ut_wo_ti = unique( tb[setdiff(1:n, i)] )
    k = length( ut_wo_ti )

    if in(ti, ut_wo_ti) # if tᵢ is not a singleton
      p = G0[b]
    end

    probs = fill(a * f(y[i], p) , k+1)
    for j in 1:k
      utj = ut_wo_ti[j]
      probs[j] = sum(utj .== tb) * f(y[i], utj)
    end

    t[b,i] = wsample( [ut_wo_ti; p] , probs )
  end # for i in 1:n
  
  # Update θ | y
  ut = unique(t[b,:])
  tb = t[b,:]
  for u in ut
    u = u[1]
    ind = find(x -> x == u, tb)
    cand = rand(Normal(u,cs))
    lg1 = sum( [lf(y[i],cand) for i in ind] ) + lg0(cand)
    lg2 = sum( [lf(y[i],u) for i in ind] ) + lg0(u)
    lg = lg1 - lg2

    if lg > log(rand())
      #acc_t += 1
      t[b,ind] = cand
    end
  end

  if b%(B/20)==0 print("\r",round(100*b/B),"%") end
  #print("\r",b)
end # for b in 1:B

#=
  include("auxGibbs.jl")
  tab = table(t[B,:])
  mt = [mean(t[:,j]) for j in 1:n]
  mns = [length(unique(t[b,:])) for b in (burn+1):B]
  tabm = table(mns)
  sort_table(tabm,true)
=#
