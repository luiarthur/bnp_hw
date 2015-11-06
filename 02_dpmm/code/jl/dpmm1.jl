include("dp.jl")
run(`mkdir -p temp`)

y1 = readdlm("../../dat/hw2.dat")

B = 100000
burn = B * .3

theta = zeros(B)
phi = ones(B)
J = 1000
G = zeros(B,J)
a = ones(B)
mu = zeros(B)
t2 = ones(B)

for b in 2:B
  
end
