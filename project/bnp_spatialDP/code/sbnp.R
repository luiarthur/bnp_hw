# https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?&ye=2004&lat=4133&submit=Submit&me=12&sitelev=&email=&step=1&p=swv_dwn&de=31&ms=1&ys=2004&plot=swv_dwn&lon=114124&ds=1
# http://www.wunderground.com/history/airport/KAAT/2000/1/1/MonthlyHistory.html?req_city=Alturas&req_statename=California&reqdb.zip=&reqdb.magic=&reqdb.wmo=&MR=1
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


gd.pts <- 70
lon.seq <- seq(-122,-120,len=gd.pts)
lat.seq <- seq(37,35,len=gd.pts)
grid <- matrix(0,gd.pts^2,3)
for (i in 1:gd.pts) {
  for (j in 1:gd.pts) {
    grid[(i-1)*gd.pts + j,1:2] <- c(lon.seq[i],lat.seq[j])
  }
}
grid[,3] <- rnorm(gd.pts^2,1:2)
grid <- grid[which(-84.5 - grid[,1]>grid[,2]),]
grid <- grid[which(-85.5 - grid[,1]<grid[,2]),]

quilt.plot(grid[,1],grid[,2],grid[,3], #lon,lat,val
           xlim=c(-122.2,-119.8), ylim=c(34.8,37.2))
map('county',add=T,lwd=2)
#abline(a=-85.5,b=-1,lwd=3)
#abline(a=-84.5,b=-1,lwd=3)


### DATA: 
# https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?&ye=2004&lat=4133&submit=Submit&me=12&sitelev=&email=&step=1&p=swv_dwn&de=31&ms=1&ys=2004&plot=swv_dwn&lon=114124&ds=1

