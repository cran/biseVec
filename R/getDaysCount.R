.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

getDaysCount <- 
function(date) {
	#Calculates how many days there are in the year of the integer date (YYMMDD)
	tmp <- 0
	days <- .C("getDaysCount", as.integer(date), d=as.integer(tmp), PACKAGE="biseVec")$d
	return(days)
}

