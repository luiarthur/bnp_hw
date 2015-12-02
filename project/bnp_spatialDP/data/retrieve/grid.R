# https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?&ye=2004&lat=4133&submit=Submit&me=12&sitelev=&email=&step=1&p=swv_dwn&de=31&ms=1&ys=2004&plot=swv_dwn&lon=114124&ds=1
# Generate Grid Points of Interest:
library(LatticeKrig) # quilt.plot
library(maps) # map()

gd.locs <- 10

lon.up <- -114#-120
lon.lo <- -126#-122
lat.up <- 44#37
lat.lo <- 32#35
lon.seq <- seq(lon.lo,lon.up,len=gd.locs)
lat.seq <- seq(lat.lo,lat.up,len=gd.locs)
int.up <- -79.5
int.lo <- -85.5

grid <- matrix(0,gd.locs^2,3)
colnames(grid) <- c("lon","lat","value")
for (i in 1:gd.locs) {
  for (j in 1:gd.locs) {
    grid[(i-1)*gd.locs + j,1:2] <- c(lon.seq[i],lat.seq[j])
  }
}
grid[,3] <- rnorm(gd.locs^2,1:2)
grid <- grid[which(int.up - grid[,1]>grid[,2]),]
grid <- grid[which(int.lo - grid[,1]<grid[,2]),]

quilt.plot(grid[,1],grid[,2],grid[,3], #lon,lat,val
           #xlim=c(lon.lo-1,lon.up-1), ylim=c(lat.lo-1,lat.up-1),
           bty='n', fg='grey',
           col= colorRampPalette(c('dark red','white','dark blue'))(100))
map('county',add=T,lwd=2,col='pink')
map('state',add=T,lwd=2,col='grey40')

write.csv(round(grid[,1:2],2),"gridlocs.dat",quote=F,row.names=F)
print(nrow(grid))
