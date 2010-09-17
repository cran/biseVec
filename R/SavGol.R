.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

SavGol <-
function(ndvi, window=7, degree=2, smoothing=10,globalthreshold=0, resvec=c(0.55,0.6,0.65,0.65,0.7)){
	days <- length(ndvi)
	res <- correctedndvi <- vector(mode="numeric",length=days);
	ndvi <- ndvi / 10000;
	
	ndvi <- ifelse(ndvi <= 0 | ndvi > 1 , NA, ndvi)
	
	#evaluate data
	if (length(which(is.na(ndvi) == FALSE)) < 10){
		return(rep(NA,7))
	}
	if (ndvi[order(ndvi, decreasing=TRUE)[10]] < 0.1) {
		return(rep(NA,7))	
	}

	ndvi <- ifelse(is.na(ndvi), -1, ndvi)

	#Savitzky-Golay Smoothing
	coef <- vector(mode="numeric", length=window);
	coef <- createCoef(window, degree);

	res <- .C("SavGol", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
		coef=as.numeric(coef),nCoef=as.integer(window) ,cndvi=as.numeric(correctedndvi),
		PACKAGE="biseVec")$cndvi

	for (i in 0:(smoothing-1)){	
		res <- .C("SavGol", rdays=as.integer(days), ndvi=as.numeric(res), 
			coef=as.numeric(coef),nCoef=as.integer(window) ,cndvi=as.numeric(correctedndvi),
			PACKAGE="biseVec")$cndvi
	}

	res <- ifelse(res <= 0, NA, res)
	
	if (length(which(is.na(res)==FALSE)) < 5){
		return(rep(NA,7))
	}

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
