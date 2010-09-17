.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

globalThresGreenup <- 
function(model, globalthreshold){
	#Green-Up Day, Global Threshold Method
	ndvi.threshold <- globalthreshold
	obsday <- NA
	model <- ifelse(is.na(model), 0, model) 
	modelhalf <- model[1:which(model == max(na.omit(model)))[1]]
	thresvec <- modelhalf-rep(ndvi.threshold, length(modelhalf))
	obsday <- which(abs(na.omit(thresvec)) == min(abs(na.omit(thresvec))))
	if (length(obsday) > 1) {
			obsday <- order(obsday, decreasing=FALSE)[1]
	}
	return(obsday)
}