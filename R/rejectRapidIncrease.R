.packageName <- "biseVec"

.First.lib <- function(lib,pkg) {
    library.dynam("biseVec", pkg, lib)
}

rejectRapidIncrease <-
function(res) {
	peaks <- c()
	threshold <- 1.5
	check <- which(is.na(res) == FALSE)
	for (i in 1:length(check)){
		if (i==1) {
			if (res[check[i]] > (threshold * mean(res[c(check[length(check)], check[i+1])]))){
				if (res[check[i]] > 0.3) {
					peaks <- c(peaks,check[i])
				}
			}
		} else {
			if (i==length(check)){
				if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[1])]))){
					if (res[check[i]] > 0.3) {
						peaks <- c(peaks,check[i])
					}
				}
			} else {
				if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[i+1])]))){
					if (res[check[i]] > 0.3) {
						peaks <- c(peaks,check[i])
					}
				}
			}
		}
	}
	res[peaks] <- NA
	return(res)
}