### Reduce Data
latlon <- read.csv("gridlocs.dat")
files <- list.files("temps")
filenum <- as.numeric( sub(".dat",files,repl="") )
L <- length(files)
ord <- order(filenum)
i <- 0
test <- read.table('temps/0.dat',skip=8,header=T) 
n <- nrow(test)

# Assuming Years are 
uyears <- unique( test$YEAR )
nyears <- length(uyears)
umonths <- 1:12
nmonths <- 12
TT <- nyears * nmonths
Y <- as.list(n)

for (f in files[ord]) {
  i <- i + 1
  dat <- read.table(paste0('temps/',f),skip=8,header=T)
  lon <- latlon[i,1]
  lat <- latlon[i,2]
  fbig <- matrix(0, TT, 5)
  colnames(fbig) <- c("year","mo","lat","lon","C")

  j <- 0
  for (yr in uyears) {
    for (mo in umonths) {
      j <- j + 1
      fbig[j,1] <- yr
      fbig[j,2] <- mo
      fbig[j,3] <- lat
      fbig[j,4] <- lon
      ind <- as.logical( (dat$YEAR == yr) * (dat$MO == mo) )
      fbig[j,5] <- mean(dat$T10M[which(ind)])
    }
  }
  Y[[i]] <- fbig
  cat("\rProgress: ", i, "/",L)
}

### Visualize
library(LatticeKrig) # quilt.plot
library(maps) # map()

view <- function(y,mo,yr) {
  M <- matrix(0,L,ncol(y[[1]]))
  colnames(M) <- colnames(y[[1]])
  i <- 0
  for (yy in y) {
    i <- i + 1
    ind <- as.logical((yy[,"year"]==yr) * (yy[,"mo"]==mo))
    M[i,] <- yy[ind,]
  }
  lon.up <- -114#-120
  lon.lo <- -124.5#-122
  lat.up <- 42#37
  lat.lo <- 32.5#35

  quilt.plot(M[,"lon"],M[,"lat"],M[,"C"], #lon,lat,val
             xlim=c(lon.lo-1,lon.up+1), ylim=c(lat.lo-1,lat.up+1),
             bty="n",fg="grey",
             col= colorRampPalette(c('dark blue','white','dark red'))(100))
  map('county',add=T,lwd=2,col='pink')
  map('state',add=T,lwd=2,col='grey40')
}

par(mfrow=c(1,3))
  view(Y,7,1985)
  view(Y,7,1994)
  view(Y,7,2004)
par(mfrow=c(1,1))

