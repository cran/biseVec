.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

biseGauss <-
function(ndvi, slidperiod=40, globalthreshold=0, asym=FALSE, resvec=c(0.55,0.6,0.65,0.65,0.7)){
	days <- length(ndvi)
	model <- res <- correctedndvi <- vector(mode="numeric",length=days)

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

	#linear interpolation (to calculate variance)
	f <- approxfun(x=1:(3*days), y=c(res,res,res))
	temp <- f((days+1):(2*days))

	##########Gaussian-Function###############
	#parameters of gaussian function
	maxpos <- order(res, decreasing=TRUE)
	maxpos <- maxpos[which(!is.na(res[maxpos]))]
	count <- 1
	mu <- maxpos[count]
	while((mu > 210)||(mu < 120)){
		count <- count + 1
		if (count > length(maxpos)){
			mu <- maxpos[1]
			break
		}
		mu <- maxpos[count]
	}
	sig <- sqrt(var(na.omit(temp)))
	base <- 0
	base <- mean(temp[1:(days/6)])

	#calculate gaussian function
	if (asym){
		res <- ifelse(is.na(res), -1, res)
		model <- .C("asymgauss", rdays=as.integer(days), ndvi=as.numeric(res),
			mustart=as.integer(mu), sigstart=as.numeric(sig), 
			rbase=as.numeric(base), model=as.numeric(model), 
			PACKAGE="biseVec")$model
	} else {
		ndvi <- res[which(is.na(res)==FALSE)]
		time <- which(is.na(res)==FALSE)
		G <- function(mu, sig, scal, base, time, days){
			erg <- ((scal / (sig * sqrt(2 * pi)))*exp(-0.5*((((time/days) - (mu/days)) / sig)^2)))+base
			return(erg)
		} 
		model.nls <- try(nls(ndvi ~ G(mu=mu, sig=sig, scal=scal, base=base, time=time, days=days), 
			start=list(sig=sig, scal=1,mu=mu), control=list(maxiter=200)),silent=TRUE)
		if (inherits(model.nls, "try-error")==FALSE){
			model <- predict(model.nls, list(time=1:days))
		} else {
			res <- ifelse(is.na(res), -1, res)
			model <- .C("gauss", rdays=as.integer(days), ndvi=as.numeric(res),
				mustart=as.integer(mu), sigstart=as.numeric(sig), 
				rbase=as.numeric(base), model=as.numeric(model), 
				PACKAGE="biseVec")$model
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
