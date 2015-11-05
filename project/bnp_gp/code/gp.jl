using Gadfly

CMAQ =  readdlm("../data/CMAQ.csv", ',', header=true)
O3 =    readdlm("../data/Ozone.csv", ',', header=true)
predL = readdlm("../data/PredLocs.csv", ',', header=true)
cmaq = CMAQ[1][:,[1,2,4]] # Lon, Lat, CMAQ_O

function matern(d, v, p, r)
  #=
    d: distance.
    v: smoothness. v increase => smoother.
    p: decay. p increase => corr (at fixed d) decreases.
  =#
  z = 2 * p  * sqrt(v) * d
  K = besselk(v,z)
  return K * z^v / ( 2^(v-1) * gamma(v) ) 
end

function GP(y, X, preds, v, init=None)
  
end
