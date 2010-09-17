.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

createCoef <-
function(window, degree){
	degree <- degree + 1
	A <- matrix(NA, window, degree)
	base <- (1:window - ceiling(window/2))
	exp <- 0:(degree-1)
	for (i in 1:window){
		for (j in 1:degree){
			A[i,j] <- base[i]^exp[j]
		}
	}
	Inv <- qr.solve(t(A) %*% A)
	E <- matrix(0, window, window)
	for (i in 1:window){
		E[i,i] <- 1
	}
	c <- vector(mode="numeric", length=window)
	for (i in 1:window){
		c[i] <- (Inv %*% (t(A) %*% E[i,]))[1]
	}
	return(c)
}