############################################################
# This script shows how to do stratified analysis using Ovary as an example
# The similar analysis can be done for other cancer types
###########################################################

library(survival)
library(ggplot2)
library(gplots)
library(ClassDiscovery)
library(ClassComparison)
library(maptools)
library(ConsensusClusterPlus)
library(superheat)
library(graphics)

## load the functions and data
setwd('path to data')
source('./CondProbMatchbyAgeGender_byYear.R')
source('./SEERdataProcess.R')
source('./timeToStable.R')
source('./CondProbGap.R')

lifeTable.female <- read.csv("./LifeTimeTable_female.csv")[,1:4]
hazard.female <- lifeTable.female[,3]
survival.female <- as.numeric(as.character(lifeTable.female[,4]))/100000
hazard.table.female <- data.frame(hazardFunc=hazard.female,survivalFunc=survival.female,age=0:100)

lifeTable.male <- read.csv("./LifeTimeTable_male.csv")[,1:4]
hazard.male <- lifeTable.male[,3]
survival.male <- as.numeric(as.character(lifeTable.male[,4]))/100000
hazard.table.male <- data.frame(hazardFunc=hazard.male,survivalFunc=survival.male,age=0:100)

subtype <- read.delim('./subtype.txt',header=F,stringsAsFactors=F)

CondProb.Diff.subtype <- lowerCI.Diff.subtype <- upperCI.Diff.subtype <- matrix(NA,nrow=4,ncol=41)
medianSurv.subtype <- sampleSize.subtype <- NULL

SEER.tmp <- read.delim('D:/SEER/SEER_1973_2013_TEXTDATA/incidence/yr1973_2013.seer9/data_tissuetype/OVARY.txt_data.txt')
SEER.tmp <- SEERdata.process(SEER.tmp)

hist.recode <- list(serous=c('8140','8440','8441','8450'),clearCell='8310',Endometrioid=c('8380','8382','8383'),Mucinous=c('8470','8471','8480','8481'))

subtype.tmp <- SEER.tmp[which(SEER.tmp$Histology_ICD_O_3 %in% hist.recode[[1]]),]
CondProb.tmp <- CondProb.Survival.control.matchedbyAgeGender.plot(subtype.tmp,'Ovary','Serous')
tumorCondProb.spline.subtype <- CondProb.tmp$CondProb.smooth 
surv.tmp <- CondProb.tmp$suvivalFunc

# diff of conditional probability in cancer and the matched normal
diff.CondProb.subtype <- tumorCondProb.spline.subtype$y-CondProb.tmp$age.gender.matched.nomralCondProb[tumorCondProb.spline.subtype$x+1]
CondProb.Diff.subtype[1,(tumorCondProb.spline.subtype$x+1)] <- diff.CondProb.subtype

lowerCI.Diff.subtype[1,] <- CondProb.tmp$CI.lower.smooth-CondProb.tmp$age.gender.matched.nomralCondProb
upperCI.Diff.subtype[1,] <- CondProb.tmp$CI.upper.smooth-CondProb.tmp$age.gender.matched.nomralCondProb
medianSurv.subtype[1] <- summary(surv.tmp)$table['median']
sampleSize.subtype[1] <- CondProb.tmp$sampleSize

pdf('C:/SEER/plots/Morta-ovary-subtype.pdf')
par(mar=c(5.1,6.1,4.1,2.1))
plot(0:40,CondProb.Diff.subtype[1,],main='SEER ovarain cancer data',
     xlab='Follow up Time (Years)',ylab='Annualized Mortality Gap',type='l',ylim=c(-0.1,0.5),
     col='magenta',lwd=2)
lines(0:40,lowerCI.Diff.subtype[1,],lty=2,col='magenta')  
lines(0:40,upperCI.Diff.subtype[1,],lty=2,col='magenta')  

colSet <- c('magenta','orangered4','cyan3','midnightblue')
for (i in 2:length(hist.recode)){
  subtype.tmp <- SEER.tmp[which(SEER.tmp$Histology_ICD_O_3 %in% hist.recode[[i]]),]
  CondProb.tmp <- CondProb.Survival.control.matchedbyAgeGender.plot(subtype.tmp,'Ovary','Serous')
  tumorCondProb.spline.subtype <- CondProb.tmp$CondProb.smooth 
  surv.tmp <- CondProb.tmp$suvivalFunc
  
  # diff of conditional probability in cancer and the matched normal
  diff.CondProb.subtype <- tumorCondProb.spline.subtype$y-CondProb.tmp$age.gender.matched.nomralCondProb[tumorCondProb.spline.subtype$x+1]
  CondProb.Diff.subtype[i,(tumorCondProb.spline.subtype$x+1)] <- diff.CondProb.subtype
  
  lowerCI.Diff.subtype[i,] <- CondProb.tmp$CI.lower.smooth-CondProb.tmp$age.gender.matched.nomralCondProb
  upperCI.Diff.subtype[i,] <- CondProb.tmp$CI.upper.smooth-CondProb.tmp$age.gender.matched.nomralCondProb
  medianSurv.subtype[i] <- summary(surv.tmp)$table['median']
  sampleSize.subtype[i] <- CondProb.tmp$sampleSize
  
  lines(0:40,CondProb.Diff.subtype[i,],lty=1,lwd=2,col=colSet[i])
  lines(0:40,lowerCI.Diff.subtype[i,],lty=2,col=colSet[i])  
  lines(0:40,upperCI.Diff.subtype[i,],lty=2,col=colSet[i])  
}

legend('topleft',legend=c(paste0(c('Serous','Clear Cell','Endometrioid','Mucinous'),', n=',sampleSize.subtype,'\nmedian survival=',
                                  round(medianSurv.subtype,1),' Years'),'95% CI'),col=c(colSet,'black'),lty=c(rep(1,4),2),y.intersp=1.7,bg='transparent',cex=0.8, bty="n")
dev.off()

