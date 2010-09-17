.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

fft.filter <-
function(ndvi, filter=8, globalthreshold=0, resvec=c(0.55,0.6,0.65,0.65,0.7)){
	#initialize
	days <- length(ndvi)
	
	ndvi <- ndvi / 10000;
	ndvi <- ifelse( ndvi <= 0 | ndvi > 1, NA, ndvi)

	#evaluate data
	if (length(which(is.na(ndvi) == FALSE)) < 10){
		return(rep(NA,7))
	}
	if (ndvi[order(ndvi, decreasing=TRUE)[10]] < 0.1) {
		return(rep(NA,7))	
	}
	
	f <- approxfun(x=1:(3*days), y=c(ndvi,ndvi,ndvi))
	ndvi.appr <- f((days+1):(2*days))
	##########################Fourier-Filter##############################
	
	ndvi.fft <- ndvi.appr

	#Fast Fourier Transfusion
	res.fft <- fft(ndvi.fft)

	#Filter
	res.filtered <- ifelse(abs(Re(res.fft)) < filter, 0, res.fft)

	#inverse
	res.inv <- fft(res.filtered, inverse=TRUE)
	res.inv.real <- vector(mode="numeric", length=days)
	for (i in 1:days){
		res.inv.real[i] <- (1/days)*sqrt(Re(res.inv[i])^2+Im(res.inv[i])^2)
	}

	model <- res.inv.real

	###########################Green-Up Day###############################

	sosvec <- vector(mode="integer", length=7)	

	#Green-Up Day, Local Threshold Method
	sosvec[1] <- localThresGreenup(model, resvec[1])
	
	#Green-Up Day, Local Threshold Method
	sosvec[2] <- localThresGreenup(model, resvec[2])
	
	#Green-Up Day, Local Threshold Method
	sosvec[3] <- localThresGreenup(model, resvec[3])

	#Green-Up Day, Global Threshold Method 0.65
	sosvec[4] <- globalThresGreenup(model, resvec[4])

	#Green-Up Day, Global Threshold Method 0.7
	sosvec[5] <- globalThresGreenup(model, resvec[5])

	#Minimum
	sosvec[6] <- min(na.omit(model))

	#Maximum
	sosvec[7] <- max(na.omit(model))

	return(sosvec)
}
