.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

biseDSig <-
function(ndvi, slidperiod=40, mode=0, resvec=c(0.55,0.6,0.65,0.65,0.7)){
	#initialize
	days <- length(ndvi)
	res <- correctedndvi <- vector(mode="numeric",length=days);
	
	ndvi <- ndvi / 10000;
	ndvi <- ifelse( ndvi <= 0 | ndvi > 1, NA, ndvi)

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

	#evaluate res
	if (length(which(is.na(res) == FALSE)) < 5){
		return(rep(NA,7))
	}
	if (res[order(res, decreasing=TRUE)[5]] < 0.1){
		return(rep(NA,7))
	}

	res <- rejectRapidIncrease(res)
	
	ndvi <- res[which(is.na(res)==FALSE)]	
	time <- which(is.na(res)==FALSE)

	maxind <- order(res, decreasing=TRUE)[1]
	#calculate double sigmoid function
	pos1 <- 0.75*maxind
	pos2 <- 1.25*maxind
	width1 <- maxind
	width2 <- days-maxind

	count <- 1
	pos1add <- -1
	pos2add <- 1
	repeat {
		model.nls <- try(nls(ndvi ~ 0.5*(tanh((time-pos1)/width1)-tanh((time-pos2)/width2)), 
			start=list(pos1=pos1, pos2=pos2, width1=width1, width2=width2),
			control=list(maxiter=200)), silent=TRUE)
		if (inherits(model.nls, "try-error")==FALSE){
			model <- predict(model.nls, list(time=1:days))
			break
		} else {
			count <- count+1
			if (count > 100){
				model <- rep(NA, days)
				break
			}
			pos1 <- pos1 + pos1add
			if (pos1 < 0){
				pos1 <- 0.75*maxind + 1
				pos1add <- 1
			}
			pos2 <- pos2+pos2add
			if (pos2 > days){
				pos2 <- 1.25*maxind - 1
				pos2add <- -1
			}
		}
	}


	if (length(which(is.na(model)==FALSE)) < 10){
		return(rep(NA,7))
	}
	
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
