latlon <- read.csv("gridlocs.dat")
files <- list.files("temps")
filenum <- as.numeric( sub(".dat",files,repl="") )
L <- length(files)
ord <- order(filenum)
i <- 0
test <- read.table(paste0('temps/',f),skip=8,header=T) 
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
      fbig[j,5] <- mean(dat$swv_dwn[which(ind)])
    }
  }
  Y[[i]] <- fbig
  cat("\rProgress: ", i, "/",L)
}
