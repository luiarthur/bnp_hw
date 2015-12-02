# https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?&ye=2004&lat=4133&submit=Submit&me=12&sitelev=&email=&step=1&p=swv_dwn&de=31&ms=1&ys=2004&plot=swv_dwn&lon=114124&ds=1
# Generate Grid Points of Interest:
library(LatticeKrig) # quilt.plot
library(maps) # map()

set.seed(1)
gd.locs <- 100

lon.up <- -114#-120
lon.lo <- -124.5#-122
lat.up <- 42#37
lat.lo <- 32.5#35
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

ind <- sample(1:nrow(grid),100)
#grid <- grid[ind,]

quilt.plot(grid[ind,1],grid[ind,2],grid[ind,3], #lon,lat,val
           xlim=c(lon.lo-1,lon.up+1), ylim=c(lat.lo-1,lat.up+1),
           bty='n', fg='grey',cex=7,
           col= colorRampPalette(c('dark red','white','dark blue'))(100))
map('county',add=T,lwd=2,col='pink')
map('state',add=T,lwd=2,col='grey40')

write.csv(round(grid[ind,1:2],2),"gridlocs.dat",quote=F,row.names=F)
write.csv(round(grid[-ind,1:2],2),"predlocs.dat",quote=F,row.names=F)
print(nrow(grid[ind,]))
