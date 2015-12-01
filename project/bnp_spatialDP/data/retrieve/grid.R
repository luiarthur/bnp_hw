# https://eosweb.larc.nasa.gov/cgi-bin/sse/daily.cgi?&ye=2004&lat=4133&submit=Submit&me=12&sitelev=&email=&step=1&p=swv_dwn&de=31&ms=1&ys=2004&plot=swv_dwn&lon=114124&ds=1
# Generate Grid Points of Interest:
library(LatticeKrig) # quilt.plot
library(maps) # map()

gd.locs <- 70
lon.seq <- seq(-122,-120,len=gd.locs)
lat.seq <- seq(37,35,len=gd.locs)
grid <- matrix(0,gd.locs^2,3)
colnames(grid) <- c("lon","lat","value")
for (i in 1:gd.locs) {
  for (j in 1:gd.locs) {
    grid[(i-1)*gd.locs + j,1:2] <- c(lon.seq[i],lat.seq[j])
  }
}
grid[,3] <- rnorm(gd.locs^2,1:2)
grid <- grid[which(-84.5 - grid[,1]>grid[,2]),]
grid <- grid[which(-85.5 - grid[,1]<grid[,2]),]
print(nrow(grid))

quilt.plot(grid[,1],grid[,2],grid[,3], #lon,lat,val
           xlim=c(-122.2,-119.8), ylim=c(34.8,37.2),bty='n',fg='grey',
           col= colorRampPalette(c('dark red','white','dark blue'))(100))
map('county',add=T,lwd=2,col='grey40')

write.csv(round(grid[,1:2],2),"gridlocs.dat",quote=F,row.names=F)
###

