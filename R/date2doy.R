.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

date2doy <- 
function(date) {
	#Calculates the Julian Day out of the integer date (YYMMDD)
	tmp <- 0
	doy <- .C("date2doy", as.integer(date), doy=as.integer(tmp), PACKAGE="biseVec")$doy
	return(doy)
}