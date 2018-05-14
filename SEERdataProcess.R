SEERdata.process <- function(SEER.tmp,...){
  remove.secondary <- which(SEER.tmp$Sequence_Number_Central > 1)
  if (length(remove.secondary) !=0){
    SEER.tmp <- SEER.tmp[-remove.secondary,]}
  removeIndex.surv1 <- grep(9999,SEER.tmp$Survival_months)
  if (length(removeIndex.surv1!=0)){
    SEER.tmp <- SEER.tmp[-removeIndex.surv1,]}
  removeIndex.surv2 <- which(SEER.tmp$Survival_months==0)
  if (length(removeIndex.surv2!=0)){
    SEER.tmp <- SEER.tmp[-removeIndex.surv2,]}
  removeIndex.age <- which(SEER.tmp$Age_at_diagnosis==999)
  if (length(removeIndex.age)!=0){
    SEER.tmp <- SEER.tmp[-removeIndex.age,]}
  over100Index <- which(SEER.tmp$Age_at_diagnosis>=100)
  if (length(over100Index)!=0){
  SEER.tmp <- SEER.tmp[-over100Index,]}

  SEER.tmp$Surv_Years <- SEER.tmp$Survival_months/12
  SEER.tmp$Vital_Status1 <- 0
  SEER.tmp$Vital_Status2 <- "Alive"
  SEER.tmp$Vital_Status1[which(SEER.tmp$Vital_Status_recode==4)] <- 1
  SEER.tmp$Vital_Status2[which(SEER.tmp$Vital_Status_recode==4)] <- "Dead"
  
#  age.at.endpoint <- floor(SEER.tmp$Age_at_diagnosis+SEER.tmp$Surv_Years)
  
#  over100Index <- which(age.at.endpoint>=100)
#  if (length(over100Index)!=0){
#    SEER.tmp <- SEER.tmp[-over100Index,]}

  return(SEER.tmp)  
}