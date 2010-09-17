.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

biseLinIP <-
function(ndvi, slidperiod=40, mode=0, resvec=c(0.55,0.6,0.65,0.65,0.7)){
	days <- length(ndvi)
	res <- vector(mode="numeric",length=days);
	correctedndvi <- vector(mode="numeric", length=days);
	
	ndvi <- ndvi / 10000;
	ndvi <- ifelse( ndvi <= 0 | ndvi > 1 ,NA, ndvi)

	#evaluate data
	if (length(which(is.na(ndvi) == FALSE)) < 10){
		return(rep(NA,7))
	}
	if (ndvi[order(ndvi, decreasing=TRUE)[10]] < 0.1) {
		return(rep(NA,7))	
	}
	
	ndvi <- ifelse( is.na(ndvi), 0, ndvi)

	#best slope index extraction
	res <- .C("bise", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
		rslidperiod=as.integer(slidperiod), cndvi=as.numeric(correctedndvi), 
		PACKAGE="biseVec")$cndvi
	
	res <- ifelse ( res <= 0, NA, res)
	
	if (length(which(is.na(res) == FALSE)) < 5){
		return(rep(NA,7))
	}

	res <- rejectRapidIncrease(res)

	#linear interpolation

	f <- approxfun(x=1:(3*days), y=c(res,res,res))
	res <- f((days+1):(2*days))

	###########################Green-Up Day###############################
	
	model <- res

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
