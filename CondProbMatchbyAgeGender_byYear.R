#################################################################################
# This function match the hazard function in general population by age and gender
# written by Yang Zhao
#################################################################################

CondProb.Survival.control.matchedbyAgeGender.plot <- function(data,tissueType,subtype=NULL,
                                                            male.hazard.table=hazard.table.female,
                                                            female.hazard.table=hazard.table.male,...){

  # survival function
  # sample size
  sampleSize.tmp <- nrow(data)
  # hazard function
  
  # calculate age at the endpoint.
  age.at.endpoint <- floor(data$Age_at_diagnosis+data$Surv_Years)
  data$Vital_Status1[which(age.at.endpoint >= 100)] <- 0
  data$Surv_Years[which(age.at.endpoint >= 100)] <- 99-data$Age_at_diagnosis[which(age.at.endpoint >= 100)]
  age.at.endpoint <- floor(data$Age_at_diagnosis+data$Surv_Years)
  age.all <- unlist(lapply(data$Age_at_diagnosis,function(x)(x+seq(0,40,1))))
  
#  if (length(which(age.all >= 100))!=0){
#  age.all[which(age.all >= 100)] <- NA}
  
  age.all <- matrix(age.all,byrow=T,ncol=41)
  for (i in 1:nrow(age.all)){
    age.all[i,which(age.all[i,]>age.at.endpoint[i])] <- NA
#	if (length(which(age.all[i,]>=100))!=0){
#	     age.all[i,which(age.all[i,] >= 100)] <- NA}
  }
  colnames(age.all) <- paste0('age.year',0:40)
  

  survival.tmp <- survfit(Surv(Surv_Years,Vital_Status1)~ 1,data=data)
  surv.estimator <- survival.tmp$surv
  surv.time <- survival.tmp$time
  tumor.CondProbbyYear.index <- unique(unlist(lapply(0:round(max(surv.time)),function(x)which.min(abs(surv.time - x)))))
  surv.time.year <- round(surv.time[tumor.CondProbbyYear.index])
  surv.time.oneYearlater <- (surv.time.year+1)
  surv.time.oneYearlater.Index <- which(surv.time.year %in% surv.time.oneYearlater)
  surv.time.Index <- surv.time.oneYearlater.Index-1

  followupTime <- surv.time.year[surv.time.Index]
  prob.death.year.bySurv <- (surv.estimator[tumor.CondProbbyYear.index][surv.time.Index]-
                           surv.estimator[tumor.CondProbbyYear.index][surv.time.oneYearlater.Index])/surv.estimator[tumor.CondProbbyYear.index][surv.time.Index]
  CondProb.smooth.spline <- smooth.spline(followupTime,prob.death.year.bySurv,spar=0.30,nknots=min(31,length(prob.death.year.bySurv)))
  
  prob.death.year.bySurv.bootstrap <- matrix(NA,ncol=(max(round(surv.time))+1),nrow=1000)
  set.seed(123)
  for (n.bootstrap in 1:1000){
      sample.bootstrap <- sample(1:nrow(data),nrow(data),replace = TRUE)
	  data.bootstrap <- data[sample.bootstrap,]
      survival.bootstrap <- survfit(Surv(Surv_Years,Vital_Status1)~ 1,data=data.bootstrap)
      surv.estimator.bootstrap <- survival.bootstrap$surv
      surv.time.bootstrap <- survival.bootstrap$time
      tumor.CondProbbyYear.index.bootstrap <- unique(unlist(lapply(0:round(max(surv.time.bootstrap)),function(x)which.min(abs(surv.time.bootstrap - x)))))
      surv.time.year.bootstrap <- round(surv.time.bootstrap[tumor.CondProbbyYear.index.bootstrap])
      surv.time.oneYearlater.bootstrap <- (surv.time.year.bootstrap+1)
      surv.time.oneYearlater.Index.bootstrap <- which(surv.time.year.bootstrap %in% surv.time.oneYearlater.bootstrap)
      surv.time.Index.bootstrap <- surv.time.oneYearlater.Index.bootstrap-1

      followupTime.bootstrap <- surv.time.year.bootstrap[surv.time.Index.bootstrap]
      prob.death.year.bySurv.bootstrap[n.bootstrap,(followupTime.bootstrap+1)] <- (surv.estimator.bootstrap[tumor.CondProbbyYear.index.bootstrap][surv.time.Index.bootstrap]-
                           surv.estimator.bootstrap[tumor.CondProbbyYear.index.bootstrap][surv.time.oneYearlater.Index.bootstrap])/surv.estimator.bootstrap[tumor.CondProbbyYear.index.bootstrap][surv.time.Index.bootstrap]
    }
    CI.upper <- apply(prob.death.year.bySurv.bootstrap,2,function(x)quantile(x,0.975,na.rm=T))
    CI.lower <- apply(prob.death.year.bySurv.bootstrap,2,function(x)quantile(x,0.025,na.rm=T))
	if (length(which(is.na(CI.upper)))!=0){
     CI.upper.smooth <- smooth.spline((0:max(round(surv.time)))[-which(is.na(CI.upper))],CI.upper[-which(is.na(CI.upper))],spar=0.30,nknots=min(31,length(CI.upper[-which(is.na(CI.upper))])))$y
     CI.lower.smooth <- smooth.spline((0:max(round(surv.time)))[-which(is.na(CI.lower))],CI.lower[-which(is.na(CI.lower))],spar=0.30,nknots=min(31,length(CI.lower[-which(is.na(CI.upper))])))$y
     CI.smooth.time <- (0:max(round(surv.time)))[-which(is.na(CI.lower))]
	}
	else{
	CI.upper.smooth <- smooth.spline((0:max(round(surv.time))),CI.upper,spar=0.3,nknots=min(31,length(CI.upper)))$y
    CI.lower.smooth <- smooth.spline((0:max(round(surv.time))),CI.lower,spar=0.3,nknots=min(31,length(CI.lower)))$y
    CI.smooth.time <- 0:max(round(surv.time))
}

  gender <- data$Sex
  gender.tab <- table(gender)  
  
  hazard.general.match.year.est <- rep(NA,41)
  for (i in 1:41){
    if (sum(table(floor(age.all[,i]),gender)) <= 5){
      next
    }
    else if (length(gender.tab)==2){
      age.frac.male <- table(floor(age.all[,i]),gender)[,"1"]/sum(table(floor(age.all[,i]),gender))
      age.frac.female <- table(floor(age.all[,i]),gender)[,"2"]/sum(table(floor(age.all[,i]),gender))
    }
    else if (length(names(gender.tab))==1 && names(gender.tab)=="1"){
      age.frac.male <- table(floor(age.all[,i]),gender)[,1]/sum(table(floor(age.all[,i]),gender))
      age.frac.female <- rep(0,length(age.frac.male))
    }
    else if (length(names(gender.tab))==1 && names(gender.tab)=="2"){
      age.frac.female <- table(floor(age.all[,i]),gender)[,1]/sum(table(floor(age.all[,i]),gender))
      age.frac.male <- rep(0,length(age.frac.female))
    }
    else {
      next
    }
    names(age.frac.female) <- names(age.frac.male) <- rownames(table(floor(age.all[,i]),gender))
    age.match.male <- match(names(age.frac.male),hazard.table.male$age)
    age.match.female <- match(names(age.frac.female),hazard.table.female$age)
    hazard.male.match.tab <- data.frame(hazard.table.male[age.match.male,],Frac=age.frac.male)
    hazard.female.match.tab <- data.frame(hazard.table.female[age.match.female,],Frac=age.frac.female)
    if (nrow(hazard.male.match.tab)!=0){
      hazard.general.match.year.est[i] <- sum(hazard.male.match.tab[,'hazardFunc']*hazard.male.match.tab[,'Frac'],na.rm=T)+sum(hazard.female.match.tab[,'hazardFunc']*hazard.female.match.tab[,'Frac'],na.rm=T)
    }
  }
  z <- list(CI.smooth.time=CI.smooth.time,CI.upper.smooth=CI.upper.smooth,CI.lower.smooth=CI.lower.smooth,
             CI.upper=CI.upper,CI.lower=CI.lower,followupTime=followupTime,CondProb.Raw=prob.death.year.bySurv,CondProb.smooth=CondProb.smooth.spline,
			 suvivalFunc=survival.tmp,sampleSize=sampleSize.tmp,age.gender.matched.nomralCondProb=hazard.general.match.year.est)
  return(z)
}
