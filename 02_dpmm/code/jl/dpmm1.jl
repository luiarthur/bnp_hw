println("Loading packages...")
using Distributions, DataFrames

y = readdlm("../../dat/hw2.dat")
n = length(y)
B = 100000

# Parameters:
theta = zeros(B,n)
phi = ones(B)
a = ones(B)
mu = zeros(B)
t2 = ones(B)

# Priors for parameters:
prior_phi = InverseGamma(2,4) # variance of data
prior_a = Gamma(1,5) # shape & scale
prior_mu = Normal(0,5) # I don't know where the mean is
prior_t2 = InverseGamma(4,3)

for b in 2:B
  # Update theta
  # Update a
  # Update phi
  # Update mu
  # Update t2
end
