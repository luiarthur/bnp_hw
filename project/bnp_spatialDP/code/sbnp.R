# http://www.ncdc.noaa.gov/cdo-web/search
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
rain <- read.csv("../data/santacruzrainfall2.csv")
head(rain[,-(1:3)])

N <- nrow(rain)
yr  <- as.numeric(substr(rain$DATE,1,4))
mo  <- as.numeric(substr(rain$DATE,5,6))
stations <- as.character(rain$STATION)
ustation <- unique(stations)

stationView <- function(stat, s.yr=0, unit=T) {
  if(s.yr==0) par(mfrow=c(4,4))
  if (s.yr == 0) s.yr  <- as.numeric(substr(stat$DATE,1,4))
  s.mo  <- as.numeric(substr(stat$DATE,5,6))


  #print(unique(s.yr))
  for (yy in sort(unique(s.yr))) {
    stat.v <- stat$TPCP[s.yr==yy]
    stat.mn <- min(stat.v)
    stat.mx <- max(stat.v)
    stat.h <- NULL
    stat.h <- stat.v
    if (unit) stat.h <- (stat.v-stat.mn)/(stat.mx-stat.mn)
    #plot(s.mo[s.yr==yy], stat.v, type='p', main=yy)
    plot(s.mo[s.yr==yy], stat.h, type='p', main=paste(stat$STATION_NAME[1],"\n",yy), xlim=c(1,12))
    if (!all(1:12 %in% s.mo[s.yr==yy])) print(paste(yy,which(!(1:12 %in% s.mo[s.yr==yy]))))
  }
  par(mfrow=c(1,1))
}

# data1: n=4, T=16; 6, 16, 26, 29
stationView(rain[which(stations==ustation[6]),])
stationView(rain[which(stations==ustation[16]),])
stationView(rain[which(stations==ustation[26]),])
stationView(rain[which(stations==ustation[29]),])

#############################################################
ss <- c(6,16,26,29)
n <- 4
TT <- nrow(rain[which(stations==ustation[29]),])
y <- matrix(0,TT*12,n)
lat <- rep(0,4)
lon <- rep(0,4)
for (i in 1:n) {
    ind <- which(stations==ustation[ss[i]])
    y[,i] <- rain$TPCP[ind]
    lat[i] <- as.numeric( as.character(rain$LAT[which(stations==ustation[ss[i]])])[1] )
    lon[i] <- as.numeric( as.character(rain$LON[which(stations==ustation[ss[i]])])[1] )
}

library(LatticeKrig)
library(geoR)
library(maps)


quilt.plot(lon,lat,c(1,2,3,4),cex=5)
map('county',add=T)

