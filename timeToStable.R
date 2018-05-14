timeToStability <- function(tumor.temp=tumor.temp,normal.tmp=normal.tmp,...){
  tumor.HazardbyYear.index <- unlist(lapply(0:40,function(x)which.min(abs(tumor.temp$time - x))))
  hazardCI.delta <- tumor.temp$lower.ci[tumor.HazardbyYear.index] - normal.tmp
  minDelta <- min(hazardCI.delta[0:20])
  hazardCI.delta.diff <- diff(hazardCI.delta)
  year.stable <- min(which(abs(hazardCI.delta.diff) < 0.003 & hazardCI.delta[-1] < (minDelta+0.02))-1)
  if (length(year.stable)==0){
    year.stable <- min(which(abs(hazardCI.delta.diff) < 0.003) -1)
  }
  return(year.stable)
}