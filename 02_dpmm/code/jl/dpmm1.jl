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

# Priors for parameters:
r_phi(a,b) = InverseGamma(2a,4) # variance of data
r_alpha(a,b) = Gamma(1,5) # shape & scale
r_eta(a,b) = Beta()
r_mu(a,b) = Normal(0,5) # I don't know where the mean is
r_t2(a,b) = InverseGamma(4,3)

for b in 2:B
  # Update theta
  # Update alpha
  # Update eta
  # Update phi
  # Update mu
  # Update t2
end
