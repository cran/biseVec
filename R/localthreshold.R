.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

localThresGreenup <- 
function(model, localthreshold=0.55){
		#Green-Up Day, Local Threshold Method
		minmodel <- min(na.omit(model))
		maxmodel <- max(na.omit(model))
		ndvi.threshold <- ((maxmodel - minmodel) * localthreshold) + minmodel
		if (is.na(ndvi.threshold)){
				return(NA)
		}
		obsday <- NA
		model <- ifelse(is.na(model), 0, model)
		modelhalf <- model[1:which(model == maxmodel)[1]]
		thresvec <- modelhalf-rep(ndvi.threshold, length(modelhalf))
		obsday <-  which(abs(na.omit(thresvec)) == min(abs(na.omit(thresvec))))
		if (length(obsday) > 1) {
			obsday <- order(obsday, decreasing=FALSE)[1]
		}
		return(obsday)
}
