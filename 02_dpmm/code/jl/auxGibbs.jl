using Distributions, DataFrames

function table(x) 
  dict = Dict(x[1] => 1)
  for u in x
    if haskey(dict, u)
      dict[u] += 1
    else
      dict[u] = 1
    end
  end
  return dict
end

y = [rand(Normal(3,.1),30);
     rand(Normal(10,.1),16); 
     rand(Normal(20,.1),26)]
#=
  F(theta) = N(theta,.1²)
  G₀ = N(0,1)
=#

n = size(y)[1]
B = 1000
t = ones(B,n)
a = 1
p = 0
G0(;n=1) = rand(Normal(0,1),n)
f(x, theta) = pdf( Normal(theta,1), x )
lf(x, theta) = logpdf( Normal(theta,1), x )
lp(x) = logpdf( Normal(0,3), x )
acc_t = 0
cs = 2

for b in 2:B
  t[b,:] = t[b-1,:]

  # Update θ_i | θ_{-i}
  for i in 1:n
    tb = t[b,:]
    ti = tb[i]
    ut_wo_ti = unique( tb[setdiff(1:n,i)] )
    k = length( ut_wo_ti )

    if in(ti, ut_wo_ti) # if tᵢ is not a singleton
      p = G0(1)[1]
    end

    probs = fill(a * f(y[i], p) / (n-1+a), k+1)
    for j in 1:k
      utj = ut_wo_ti[j]
      probs[j] = sum(utj .== tb) * f(y[i], utj) / (n-1+a)
    end

    t[b,i] = wsample([ut_wo_ti; p], probs)
  end # for i in 1:n
  
  # Update θ | y
  ut = unique(t[b,:])
  tb = t[b,:]
  for u in ut
    u = u[1]
    ind = find(x -> x == u, tb)
    cand = rand(Normal(u,cs))
    lg1 = sum( [lf(i,cand) for i in ind] ) + lp(cand)
    lg2 = sum( [lf(i,u) for i in ind] ) + lp(u)
    lg = lg1 - lg2

    if lg > log(rand())
      #acc_t += 1
      t[b,i] = cand
    end
  end

  if b%(B/10)==0 print("\r",round(100*b/B),"%") end
end # for b in 1:B

#=
  include("auxGibbs.jl")
  tab = table(t[2,:])
  mt = [mean(t[:,j]) for j in 1:n]
=#
