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
             bty="n",fg="grey",breaks=seq(0,32,len=101),main=paste0(mo,"/",yr),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100))
  map('county',add=T,lwd=2,col='pink')
  map('state',add=T,lwd=2,col='grey40')
  M
}

#par(mfrow=c(5,4),mar=c(4,4,1,1))
par.opts <- par()
par(mfrow=c(2,3),mar=c(4,4,1,1))
  for (yr in c(1985,1989,1994,1999,2001,2004)) view(Y,10,yr)
par(mfrow=c(1,1))

Yout <- array(0,c(L,TT,5))
for (i in 1:L) { 
  Yout[i,,] <- Y[[i]]
}

save(Yout,file="Y.RData")
