###########################################################################
########## ULSAM-70, observational associations ###########################

  ## upload Metabo file
dir <- ("... /Scripts")
setwd(dir)
U70_clamp_met_nonDM <- read.table("ulsam_obs_master.txt")
 
  ## Normality diagnostic plots ##
par(mfrow=c(3,1))
hist(U70_clamp_met_nonDM$z620,prob=T,breaks=20,main="U70 n=954 nonDM Clamp",col="gray")
qqnorm((U70_clamp_met_nonDM$z620),main="U70 n=954 nonDM Clamp",col="red")
qqline((U70_clamp_met_nonDM$z620))
boxplot(U70_clamp_met_nonDM$z620,horizontal=T, col="lightblue",main="U70 n=954 nonDM Clamp")

#### 1) clamp M/I 
  ## age + quality
colids <- grep("^X0",colnames(U70_clamp_met_nonDM))
attach(U70_clamp_met_nonDM)
METABO_beta_1 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,1])
METABO_se_1 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,2])
METABO_p_1 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,4])
METABO_fstat_1 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_1 <- sapply(colids, function(x) pf(
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F)) # model F-statistic extraction
Metanames <- colnames(U70_clamp_met_nonDM[colids])
P_fdr <- (p.adjust(METABO_p_1,method="fdr")) 
LM_clamp_age_results <- cbind(Metanames, METABO_beta_1,METABO_se_1,METABO_p_1,P_fdr,METABO_fstat_1,METABO_fstat_p_1)  

  ## Residual plots for empty model
attach(U70_clamp_met_nonDM)
LM_clamp_age <- lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality )
par(mfrow=c(3,1))
hist(LM_clamp_age$res,main="U70 residuals Clamp=age+quality")
plot(LM_clamp_age, which=1)
plot(LM_clamp_age, which=2)
  
  ## output table with matched metabolite names
LM_clamp_age_results <- data.frame(LM_clamp_age_results)
ids <- gsub("^X0_","",LM_clamp_age_results$Metanames)
matchlist <- data.frame(read.delim("~/Desktop/Project III_Metabolomic Fingerprint/metabo_id_match_list.txt"))
matchid <- (match(ids,matchlist[,3]))
matchlist$matchid = 1:192
matchlist <- matchlist[,c(1:3,7)]
LM_clamp_age_results$matchid <- matchid
m <- merge(matchlist,LM_clamp_age_results,by="matchid")
m2 <- with(m,m[order(m[,2]),])

  ## additional BMI adjustement
attach(U70_clamp_met_nonDM)
METABO_beta_2 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,1])
METABO_se_2 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,2])
METABO_p_2 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,4])
METABO_fstat_2 <- sapply(colids, function(x) 
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_2 <- sapply(colids, function(x) pf(
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F))
P_fdr_2 <- (p.adjust(METABO_p_2,method="fdr"))
LM_clamp_age_bmi_results <- cbind (Metanames, METABO_beta_2,METABO_se_2,P_fdr_2,METABO_p_2,METABO_fstat_2,METABO_fstat_p_2) 

  ## Residual plots for empty model
attach(U70_clamp_met_nonDM)
LM_clamp_age_bmi <- lm(scale(Clamp) ~ Age + Storage_time + Thawed + Sample_type + Quality +BMI )
par(mfrow=c(3,1))
hist(LM_clamp_age_bmi$res,main="U70 residuals Clamp=bmi+age+quality")
plot(LM_clamp_age_bmi, which=1)
plot(LM_clamp_age_bmi, which=2)

LM_clamp_age_bmi_results <- data.frame(LM_clamp_age_bmi_results)
ids <- gsub("^X0_","",LM_clamp_age_bmi_results$Metanames)
LM_clamp_age_bmi_results$matchid <- matchid
ma <- merge(matchlist,LM_clamp_age_bmi_results,by="matchid")
ma2 <- with(ma,ma[order(ma[,2]),])

  ## Diagnostics for outliers, leverage 
  ## Cook's distance
LM_baseline <- LM_clamp_age
par(mfrow=c(2,1))
plot(LM_baseline, which=(4),main="nonDM - Cook's D (rawClamp)",cex.main=2)
  ## Studentised Residuals 
library(car)
outlierTest(LM_baseline) # studentized residual Outlier test, Bonferroni-adjusted
qqPlot(LM_baseline,simulate=T,labels=row.names(data.frame(Clamp)),main="Studentised Residuals")
  ## Influence Plot - combined stRes, Cook's, hat values (Leverage)
par(mfrow=c(1,1))
plot(hatvalues(LM_baseline),rstudent(LM_baseline),type='n',ylim=c(-3,6),
     xlab="hat values (leverage/extreme)",ylab="student. residuals (aberration from regression line)",
     main="Influence plot (circle diameter = Cook's Dist.")
cook <- sqrt(cooks.distance(LM_baseline))
points(hatvalues(LM_baseline),rstudent(LM_baseline),cex=10*cook/max(cook))
abline(v=3/25,lty=2)
abline(h=c(-2,0,2),lty=2)
Cllm <- data.frame(cbind(Clamp,Age,Storage_time,Thawed,Sample_type,Quality))
Cll <- Cllm[!is.na(Cllm$Clamp),]
identify(hatvalues(LM_baseline),rstudent(LM_baseline),row.names(Cll))

hist(Cll$Clamp,breaks=30,col="gray")
identify(rep(1,length(Cll$Clamp)),Cll$Clamp,labels=seq_along(Cll$Clamp))

#### 2) Insulinogenic index 

hist((U70_clamp_met_nonDM$IGI30))
  
  ## distrubition skewed - needs ln-transformation 
U70_clamp_met_nonDM$lnIGI30 <- log(U70_clamp_met_nonDM$IGI30)
  ## Normality, Boxplot ##
par(mfrow=c(3,1))
hist((U70_clamp_met_nonDM$lnIGI30),prob=T,breaks=20,main="U70 ln(IGI30) nonDM(red), diabetics (blue)",col="brown1",)
qqnorm((U70_clamp_met_nonDM$lnIGI30),main="U70 n=954 nonDM Clamp",col="red")
qqline((U70_clamp_met_nonDM$lnIGI30))
boxplot(U70_clamp_met_nonDM$lnIGI30,horizontal=T, col="lightblue")

attach(U70_clamp_met_nonDM)
METABO_beta_3 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,1])
METABO_se_3 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,2])
METABO_p_3 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,4])
METABO_fstat_3 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_3 <- sapply(colids, function(x) pf(
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F)) 
P_fdr_3 <- (p.adjust(METABO_p_3,method="fdr"))
LM_lnIGI30_age_results <- data.frame(cbind(Metanames, METABO_beta_3,METABO_se_3,METABO_p_3,P_fdr_3,METABO_fstat_3,METABO_fstat_p_3))
ids<-gsub("^X0_","",LM_lnIGI30_age_results$Metanames)
LM_lnIGI30_age_results$matchid<-matchid
mb<-merge(matchlist,LM_lnIGI30_age_results,by="matchid")
mb2<-with(mb,mb[order(mb[,2]),])

attach(U70_clamp_met_nonDM)
METABO_beta_4 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,1])
METABO_se_4 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,2])
METABO_p_4 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,4])
METABO_fstat_4 <- sapply(colids, function(x) 
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_4 <- sapply(colids, function(x) pf(
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(lnIGI30) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F)) 
P_fdr_4 <- (p.adjust(METABO_p_4,method="fdr"))
LM_lnIGI30_age_bmi_results <- data.frame(cbind (Metanames, METABO_beta_4,METABO_se_4,METABO_p_4,P_fdr_4,METABO_fstat_4,METABO_fstat_p_4))
ids <- gsub("^X0_","",LM_lnIGI30_age_bmi_results$Metanames)
LM_lnIGI30_age_bmi_results$matchid <- matchid
mc <- merge(matchlist,LM_lnIGI30_age_bmi_results,by="matchid")
mc2 <- with(mc,mc[order(mc[,2]),])

### 3) ORAL DISPOSITION INDEX 
hist(U70_clamp_met_nonDM$DI)
 
 ## DI needs log-TF
U70_clamp_met_nonDM$lnDI <- log(U70_clamp_met_nonDM$DI)

par(mfrow=c(3,1))
hist((U70_clamp_met_nonDM$lnDI),prob=T,breaks=20,main="U70 nonDM ln-DI",col="gray")
qqnorm((U70_clamp_met_nonDM$lnDI),main="U70 n=954 nonDM Clamp",col="red")
qqline((U70_clamp_met_nonDM$lnDI))
boxplot(U70_clamp_met_nonDM$lnDI,horizontal=T, col="lightblue")

attach(U70_clamp_met_nonDM)
METABO_beta_5 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,1])
METABO_se_5 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,2])
METABO_p_5 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$coef[7,4])
METABO_fstat_5 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_5 <- sapply(colids, function(x) pf(
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F)) 
P_fdr_5 <- (p.adjust(METABO_p_5,method="fdr"))
LM_lnDI_age_results <- data.frame(cbind (Metanames, METABO_beta_5,METABO_se_5,METABO_p_5,P_fdr_5,METABO_fstat_5,METABO_fstat_p_5))
ids <- gsub("^X0_","",LM_lnDI_age_results$Metanames)
LM_lnDI_age_results$matchid <- matchid
md <- merge(matchlist,LM_lnDI_age_results,by="matchid")
md2 <- with(md,md[order(md[,2]),])

attach(U70_clamp_met_nonDM)
METABO_beta_6 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,1])
METABO_se_6 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,2])
METABO_p_6 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$coef[8,4])
METABO_fstat_6 <- sapply(colids, function(x) 
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1])
METABO_fstat_p_6 <- sapply(colids, function(x) pf(
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[1],
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[2],
  summary(lm(scale(lnDI) ~ Age + Storage_time + Thawed + Sample_type + Quality + BMI + U70_clamp_met_nonDM[,x]))$fstatistic[3],
  lower.tail=F)) # model F-statistic extraction
P_fdr_6 <- (p.adjust(METABO_p_6,method="fdr"))
LM_lnDI_age_bmi_results <- data.frame(cbind(Metanames, METABO_beta_6,METABO_se_6,METABO_p_6,P_fdr_6,METABO_fstat_6,METABO_fstat_p_6))
ids <- gsub("^X0_","",LM_lnDI_age_bmi_results$Metanames)
LM_lnDI_age_bmi_results$matchid <- matchid
me <- merge(matchlist,LM_lnDI_age_bmi_results,by="matchid")
me2 <- with(me,me[order(me[,2]),])


