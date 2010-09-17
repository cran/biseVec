.packageName <- "bise"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

biseDSigC <-
function(ndvi, slidperiod=40, globalthreshold=0,localthreshold=0.55, draw=FALSE, color="blue",resvec=c(0.55,0.6,0.65,0.65,0.7)){
	#initialize
	days <- length(ndvi)
	res <- correctedndvi <- vector(mode="numeric",length=days);
	
	ndvi <- ndvi / 10000;
	ndvi <- ifelse( ndvi <= 0 | ndvi > 1, NA, ndvi)

	#evaluate data
	if (length(which(is.na(ndvi) == FALSE)) < 10){
		return(NA)
	}
	if (ndvi[order(ndvi, decreasing=TRUE)[10]] < 0.1) {
		return(NA)	
	}
	
	ndvi <- ifelse( is.na(ndvi), 0, ndvi)

	#best slope index extraction
	res <- .C("bise", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
		rslidperiod=as.integer(slidperiod), cndvi=as.numeric(correctedndvi), 
		PACKAGE="biseVec")$cndvi
	
	res <- ifelse ( res <= 0, NA, res)

	#evaluate res (not enough values)
	if (length(which(is.na(res) == FALSE)) < 5){
		return(NA)
	}
	if (res[order(res, decreasing=TRUE)[5]] < 0.1){
		return(NA)
	}

	res <- rejectRapidIncrease(res)
	
	if (draw){
		points(res, col="red", lwd=2)
	}

	#linear interpolation
	f <- approxfun(x=1:(3*days), y=c(res,res,res))
	res <- f((days+1):(2*days))

	###############################################################################
	
	doublesigmoid <- function(res) {
		#########################single sigmoid############################
		sigmoid <- function(res, maxvalue=0, turnaround=FALSE) {
			days <- length(res)
			model <- vector(mode="numeric", length=days)
	
			if (turnaround){
				res <- res[days:1]
			}
		
			#calculate coefficients F and G
			#F - Base, (G+F) Maximum
			F <- mean(na.omit(res[1:(days/4)]))
			#F <- 0.2
			if (maxvalue > 0){
				G <- maxvalue - F
			} else {
				res.order <- order(na.omit(res), decreasing=TRUE)
				meanmax <- mean(res[res.order[1:3]])
				G <- meanmax - F
			}

			model <- .C("sigmoid", rdays=as.integer(days), ndvi=as.numeric(res),
				rF=as.numeric(F), rG=as.numeric(G), model=as.numeric(model), 
				PACKAGE="biseVec")$model

			if (turnaround){
				model <- model[length(model):1]
			}
			if (length(model) != days){
				model <- rep(NA, days)
			}
			return(model)
		}
		##############################################################
		days <- length(res)
		delay <- 10
		count <- 1
		while((maxpos <- order(res, decreasing=TRUE)[count]) > 300){
			count <- count + 1
			if ((count > 100) || (is.na(maxpos))){
				return(rep(NA, days))
			}
		}
		#calculate front sigmoid
		modelfront <- sigmoid(res[1:(maxpos + delay)])
		if (length(na.omit(modelfront)) < 5){
			return(rep(NA, days))
		}
		#calculate back sigmoid
		maxvalue <- max(na.omit(modelfront))
		modelback <- sigmoid(res[(maxpos-1 + delay):days], maxvalue=maxvalue, turnaround=TRUE)
		if ((length(modelback) < 5) || is.na(modelback[3])){
			modelback <- rep(maxvalue, days-length(modelfront)+2)
		}	

		model <- vector(mode="numeric", length=days)

		#combine functions
		endfront <- length(modelfront)
		endback <- length(modelback)
		if (modelfront[endfront] != modelback[3]){
			scal <- modelback[3]
			for (i in 1:endback){
				modelback[i] <- (modelback[i] / scal) * modelfront[endfront]
			}
		}
		model <- c(modelfront, modelback[3:endback])
		
		return(model)
	}

	###############################################################################

	maxind <- order(res, decreasing=TRUE)[1]
	#calculate double sigmoid function
	model <- doublesigmoid(res)

	if (length(is.na(model)==FALSE) < 10){
		return(NA)
	}
	if (draw){
		lines(model, col=color, lwd=2)
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
