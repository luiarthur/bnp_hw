println("Loading packages...")
using Distributions, DataFrames

y = readdlm("../../dat/hw2.dat")
n = length(y)
B = 100000

# Parameters:
theta = zeros(B,n)
phi = ones(B)
alpha = ones(B)
eta = ones(B) # auxiliary variable
mu = zeros(B)
t2 = ones(B)

a_phi, b_phi = 2,4
a_alpha, b_alpha = 1,5
a_mu, b_mu = 0,5
a_t2, b_t2 = 4,3

# Posterior Generator for parameters:
# theta generator
function r_theta(theta_curr, alpha_curr, phi_curr, mu_curr, t2_curr)
  theta_new = theta_curr + 0

  # Update θ_i | θ_{-i}
  for i in 1:n
    q0 = (2*pi*(t2_curr+phi_curr))^-.5 * 
      exp(( (y[i] * t2_curr + mu_curr * phi_curr) / (t2_curr+phi_curr) )^2 - 
      (y[i]^2 * t2_curr + mu_curr^2 * phi_curr) / (2 * t2_curr * phi_curr)) 

    ut_xi = unique( theta_new[ setdiff(1:n,i) ] )
    nstar_xi = length( ut_xi )

    #n_xi = SharedArray( Float64, nstar_xi )
    #q = SharedArray( Float64, nstar_xi )
    #@parallel for j in 1:nstar_xi

    n_xi = zeros( nstar_xi )
    q = zeros( nstar_xi )
    for j in 1:nstar_xi
      n_xi[j] = sum(ut_xi[j] .== theta_new)
      q[j] = pdf( Normal(ut_xi[j], phi_curr) , y[i])
    end

    p_draw_new_theta = alpha_curr * q0
    p_draw_existing = n_xi .* q

    probs = [ p_draw_new_theta; p_draw_existing ]
    ind = wsample( 0:nstar_xi , probs )

    if ind == 0
      new_mean = (y[i] * t2_curr + mu_curr * phi_curr) / (t2_curr + phi_curr)
      new_sd = t2_curr * phi_curr / (2 * t2_curr * phi_curr)
      theta_new[i] = rand(Normal( new_mean , new_sd ))
    else
      theta_new[i] = ut_xi[ind]
    end
  end # for i in 1:n

  # Update θ | y
  # NEED TO DO THIS!!! LEFT OFF HERE!!!

  theta_new
end

# phi generator
function r_phi(theta_curr) 
  new_a = n/2 +a_phi
  new_b = sum([ (y[i] - theta_curr[i])^2 for i in 1:n]) / 2 + b_phi
  rand( InverseGamma(new_a , new_b) ) # variance of data
end

# mu generator
function r_mu(theta_curr, t2_curr)
  ut = unique(theta_curr)
  nstar = length(ut)
  denom = nstar * b_mu + t2_curr
  new_mean = (nstar * mean(theta_curr) * b_mu + a_mu * t2_curr) / denom
  new_sd = sqrt( b_mu*t2_curr / denom  )
  rand( Normal(new_mean, new_sd)  )
end

# tau^2 generator
function r_t2(theta_curr, mu_curr)
  ut = unique(theta_curr)
  nstar = length(ut)
  new_a = nstar / 2 + a_t2
  new_b = sum([ (ut[j]-mu_curr) ^2 for j in 1:nstar]) / 2 + b_t2
  rand( InverseGamma(new_a, new_b) )
end

# eta generator. Auxiliary variable.
r_eta(alpha_curr) = rand( Beta(alpha_curr+1,n) )

# alpha generator
function r_alpha(theta_curr, eta_curr) 
  nstar = length( unique(theta_curr) )
  c = a_alpha + nstar
  d = b_alpha - log(eta_curr)
  #eps = (c - 1) / (n*d + c - 1)
  ind = wsample([0,1],[c-1, n*d])
  out = 0 
  if ind == 0 
    out = Gamma(c, 1/d) # shape = c, and scale = 1/d => rate = d.
  else
    out = Gamma(c-1, 1/d)
  end
  rand(out)
end


# Gibbs. 8 seconds for 100 mcmc iterations.
println("Starting Computation...")
@time for b in 2:B
  theta[b,:] = r_theta(theta[b-1,:], alpha[b-1], phi[b-1], mu[b-1], t2[b-1])
  alpha[b] = r_alpha(theta[b,:], eta[b-1])
  eta[b] = r_eta(alpha[b])
  phi[b] = r_phi(theta[b,:]) 
  mu[b] = r_mu(theta[b,:], t2[b-1])
  t2[b] = r_t2(theta[b,:], mu[b])

  #if b%(B/100)==0 print("\r",round(100*b/B),"%") end
  print("\r Progress: ", b)
end

#=
  include("dpmm1.jl")
  writedlm("temp/out_alpha.dat", alpha)
  writedlm("temp/out_eta.dat", eta)
  writedlm("temp/out_phi.dat", phi)
  writedlm("temp/out_mu.dat", mu)
  writedlm("temp/out_t2.dat", t2)
  writedlm("temp/out_theta.dat", theta)
=#
