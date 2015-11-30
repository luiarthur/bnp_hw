# http://www.ncdc.noaa.gov/cdo-web/search
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
rain <- read.csv("../data/santacruzrainfall.csv")
head(rain[,-(1:3)])

n <- nrow(rain)
yr  <- as.numeric(substr(rain$DATE,1,4))
mo  <- as.numeric(substr(rain$DATE,5,6))
stations <- as.character(rain$STATION)
ustation <- unique(stations)

stationView <- function(stat, s.yr=0) {
  if(s.yr==0) par(mfrow=c(4,4))
  if (s.yr == 0) s.yr  <- as.numeric(substr(stat$DATE,1,4))
  s.mo  <- as.numeric(substr(stat$DATE,5,6))


  for (yy in sort(unique(s.yr))) {
    stat.v <- stat$TPCP[s.yr==yy]
    stat.mn <- min(stat.v)
    stat.mx <- max(stat.v)
    #plot(s.mo[s.yr==yy], stat.v, type='p', main=yy)
    plot(s.mo[s.yr==yy], (stat.v-stat.mn)/(stat.mx-stat.mn), type='p', main=paste(stat$STATION_NAME[1],"\n",yy), xlim=c(1,12))
  }
  par(mfrow=c(1,1))
}

stationView(rain[which(stations==ustation[29]),]) # 6, 16, 26, 29
stationView(rain[which(stations==ustation[16]),]) # 6, 16, 26, 29
stationView(rain[which(stations==ustation[26]),]) # 6, 16, 26, 29
stationView(rain[which(stations==ustation[6]),]) # 6, 16, 26, 29

# n=4, T=16
