############################################################################
###### MR analysis for Insulin Resistance - omit ULSAM #####################
############################################################################

dir <- (" ... /Scripts")
setwd(dir)

library(plyr)
library(metafor)
library(gtx)

  ## starting point: 52 metabolites associated insulin response
  ## master file uploads Pivus and TwinGene
P70 <- read.table("Pivus70_IR.txt")
TWGE <- read.table("TWGE_IR.txt")

  ## name age and sex similarly
TWGE$AGE <- TWGE$age
TWGE$SEX <-  - (TWGE$sex -2) # 1 is female as in PIVUS (0 is male)
P70$AGE <- P70$age70.x
P70$SEX <- P70$kn

  ## merge
COMB <- merge(TWGE,P70,all.x=T,all.y=T)
COMB$Cohort_1 <- 0
COMB$Cohort_1[which(is.na(COMB$piusnr))] <- 1 # 1 is TWGE

# write.table(COMB,"MR_P_U_for_IR.txt")

metabonames <- c("MAG_18_2","MAG_16_1","MAG_18_1","Piperine","Caffeine","LPE_18_2","DOCA_Glycine_Conj","LPE_18_1","MAG_16_0",
                 "Corticosterone","LPC_18_3","LPC_20_0","gamma_Glu_Leu","MAG_14_0","LPC_18_2_0_0","Propanolol","LPC_0_0_18_2",
                 "DAG_34_1","Arachidonic_acid_ethyl_ester","PC_36_4a","Eicosatrienoic_acid","LPE_18_1b","PC_16_2","gamma_Tocopherol",
                 "Cinnamic_acid","LPC_18_1a","LPC_18_1b","Paraxanthine_Theophylline","L_Tyrosine","Hippuric_acid","linoleoyl_2_stearoyl_sn_glycerol",
                 "Chenodeoxycholic_acid_Glyc_Conj","Decanoyl_L_carnitine","LPE_20_4","LPE_18_1_P","Glycocholic_acid","L_Acetylcarnitine","Betaine",
                 "Trihydroxy_5b_cholanoic_acid","PC_38_7","PC_38_3","Oleic_acid","LPE_16_0a","Choline","LPC_20_3a","Treprostinil","Hexose","Docosapentaenoic_acid",
                 "Bilirubin","Bilirubin_IXa","Palmitoleic_acid","Limoleic_acid")

#### MR analysis: GRS x metabolite associations, adj. for age, sex, cohort
  ## scaled metabolites, mean = 0, SD = 1

### as requested by reviewers: control for population stratification
### adjustment for the first 3 genetic principle components 

#### additional adjustment for genetic principle components PC1, PC2, PC3 #######

pivus_PC <- read.table("~/pivus.+.qced.pca.evec.txt", quote="\"")
pivus_PC$lpnr <- as.numeric(pivus_PC[,1])
twge_PC <- read.table("~/twge.+.qced.pca.evec.txt", quote="\"")
split_1 <- str_split_fixed(as.character(twge_PC[,1]), ":", 2)
twge_PC$genetic_ID <- (split_1[,1])

COMB$PC1<-c(); COMB$PC2<-c(); COMB$PC3<-c()
for (i in 1:nrow(COMB)){
  if (!is.na(COMB$lpnr[i])) COMB$PC1[i] <- pivus_PC[match(COMB$lpnr[i],pivus_PC$lpnr),2]
  else COMB$PC1[i] <- twge_PC[match(COMB$genetic_ID[i],twge_PC$genetic_ID),2]
  if (!is.na(COMB$lpnr[i])) COMB$PC2[i] <- pivus_PC[match(COMB$lpnr[i],pivus_PC$lpnr),3]
  else COMB$PC2[i] <- twge_PC[match(COMB$genetic_ID[i],twge_PC$genetic_ID),3]
  if (!is.na(COMB$lpnr[i])) COMB$PC3[i] <- pivus_PC[match(COMB$lpnr[i],pivus_PC$lpnr),4]
  else COMB$PC3[i] <- twge_PC[match(COMB$genetic_ID[i],twge_PC$genetic_ID),4]
}
COMB$PC1 <- ifelse(is.na(COMB$PC1),mean(COMB$PC1,na.rm=T),COMB$PC1) # few NA's imputed with mean  
COMB$PC2 <- ifelse(is.na(COMB$PC2),mean(COMB$PC2,na.rm=T),COMB$PC2)
COMB$PC3 <- ifelse(is.na(COMB$PC3),mean(COMB$PC3,na.rm=T),COMB$PC3)


CNT2_IR_Beta <- sapply(which(names(COMB) %in% metabonames), function(x) 
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX  + COMB$Cohort_1 + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$coef[5,1])
CNT2_IR_SE <- sapply(which(names(COMB) %in% metabonames), function(x) 
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1 + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$coef[5,2])
CNT2_IR_P <- sapply(which(names(COMB) %in% metabonames), function(x) 
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1 + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$coef[5,4])
CNT2_IR_Fstat <- sapply(which(names(COMB) %in% metabonames), function(x) 
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1 + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$fstatistic[1])
CNT2_IR_Fstat_P <- sapply(which(names(COMB) %in% metabonames), function(x) pf(
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1  + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$fstatistic[1],
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1  + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$fstatistic[2],
  summary(lm(scale(COMB[,x]) ~ COMB$AGE + COMB$SEX + COMB$Cohort_1 + COMB$CNT2 + COMB$PC1 + COMB$PC2 + COMB$PC3))$fstatistic[3],
  lower.tail=F)) # model F-statistic extraction
CNT2_IR_results <- cbind (names(COMB)[which(names(COMB) %in% metabonames)], CNT2_IR_Beta,CNT2_IR_SE,CNT2_IR_P,CNT2_IR_Fstat,CNT2_IR_Fstat_P)
temp0 <- (CNT2_IR_results[order(as.numeric(CNT2_IR_results[,4])),])

  ## IV analysis 
  ## GRS x risk factor association reported by Scott et al.
  ## the negative insulin sensitivity (clamp) association is *(-1) to indicate higher GRS == worse IR
SCOTTclamp_beta <- 0.03 
SCOTTclamp_se <- (-0.01 - -0.04)/3.92

IR_MR_beta <- sapply(1:52, function(x) CNT2_IR_Beta[x] / SCOTTclamp_beta)
IR_MR_se <- sapply(1:52, function(x) abs(IR_MR_beta[x])*sqrt( ((CNT2_IR_SE[x]/CNT2_IR_Beta[x])^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) ))
IR_MR_z <- sapply(1:52, function(x) IR_MR_beta[x] / IR_MR_se[x])
IR_MR_p <- sapply(1:52, function(x) 2*pnorm(-abs(IR_MR_z[x])))
IR_MR_results <- cbind(names(COMB)[which(names(COMB) %in% metabonames)],IR_MR_beta,IR_MR_se,IR_MR_z,IR_MR_p)
temp1 <- IR_MR_results

colnames(temp0)[1]<- "metabo"; colnames(temp1)[1]<- "metabo"
temp3 <- merge(temp0,temp1,by="metabo")
temp3 <- (temp3[order(as.numeric(temp3$IR_MR_p)),])

# write.table(temp3,"IR_MR_U70_P70_TWGE_results.txt")

##################################################################################
  ## MR replication in KORA/TwinsUK, CHARGE, Finnish study

ir_snps <- data.frame(read.delim("scott_ir_snps.txt")) # Scott et al. Suppl Table
ir_snps <- ir_snps[order(as.character(ir_snps$rs_id)),] # alph order rs_id

#### KORA/TwinsUK: available are
# M33447 - palmitoleate, M01359 - oleate, M01299 - tyrosine, M15753 - Hippurate,
# M27447 - MAG18:2 (2-monolinolein), M21184 - MAG(18:1) (1-monoolein), M33420- gamma-tocopherol

palmit_kora <- read.delim("M33447.metal.pos.txt")
oleate_kora <- read.delim("M01359.metal.pos.txt")
tyrosine_kora <- read.delim("M01299.metal.pos.txt")
hippuric_kora <- read.delim("M15753.metal.pos.txt")
mag182_kora <- read.delim("M27447.metal.pos.txt")
mag181_kora <- read.delim("M21184.metal.pos.txt")
gtoco_kora <- read.delim("M33420.metal.pos.txt")

  ## KORA: metaboltes on log10-scale - to compare studies convert to SD-unit
  ## SDs per metabolite reported in Suppl Table separately for TwinsUK and KORA
  # palimtoleate      0.208 (n=6010)      0.208 (n=1766)
  # oleate            0.140 (n=6003)      0.111 (n=1765)
  # hippurate         0.280 (n=6048)      0.294 (n=1758)
  # tyrosine          0.093 (n=6046)      0.077 (n=1761)
  # mag182 (2-monol)  0.246 (n=1039)      0.221 (n=1758)
  # mag181            0.262 (n=4027)      0.250 (n=1690)            
  # gammatoco         0.244 (n=6226)      0.218 (n=977)

  ## Var and SD as sample-size weighted average 
k_palmit_sd <- sqrt(((0.208^2)*6010 + (0.208^2)*1766)/(6010+1766))
k_oleate_sd <- sqrt(((0.140^2)*6003 + (0.111^2)*1765)/(6003+1765))
k_hippuric_sd <- sqrt(((0.280^2)*6048 + (0.294^2)*1758)/(6048+1758))
k_tyrosine_sd <- sqrt(((0.093^2)*6046 + (0.077^2)*1761)/(6046+1761))
k_mag182_sd <- sqrt(((0.246^2)*1039 + (0.221^2)*1758)/(1039+1758))
k_mag181_sd <- sqrt(((0.262^2)*4027 + (0.250^2)*1690)/(4027+1690))
k_gtoco_sd <- sqrt(((0.244^2)*6226 + (0.218^2)*977)/(6226+977))

  ## select IR snps from metabolite GWAS results and order by rs_id for later merging
oleate_kora_2 <- oleate_kora[which(as.character(oleate_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
oleate_kora_2 <- oleate_kora_2[order(as.character(oleate_kora_2$MarkerName)),] # order alph by rs_id
palmit_kora_2 <- palmit_kora[which(as.character(palmit_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
palmit_kora_2 <- palmit_kora_2[order(as.character(palmit_kora_2$MarkerName)),] # order alph by rs_id
tyrosine_kora_2 <- tyrosine_kora[which(as.character(tyrosine_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
tyrosine_kora_2 <- tyrosine_kora_2[order(as.character(tyrosine_kora_2$MarkerName)),]
hippuric_kora_2 <- hippuric_kora[which(as.character(hippuric_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
hippuric_kora_2 <- hippuric_kora_2[order(as.character(hippuric_kora_2$MarkerName)),] # order alph by rs_id
mag182_kora_2 <- mag182_kora[which(as.character(mag182_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
mag182_kora_2 <- mag182_kora_2[which(as.character(mag182_kora_2$MarkerName) %in% as.character(ir_snps$rs_id)),]
mag181_kora_2 <- mag181_kora[which(as.character(mag181_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
mag181_kora_2 <- mag181_kora_2[which(as.character(mag181_kora_2$MarkerName) %in% as.character(ir_snps$rs_id)),]
gtoco_kora_2 <- gtoco_kora[which(as.character(gtoco_kora$MarkerName) %in% as.character(ir_snps$rs_id)),]
gtoco_kora_2 <- gtoco_kora_2[which(as.character(gtoco_kora_2$MarkerName) %in% as.character(ir_snps$rs_id)),]

  ## align KORA/TwinsUK-coded with IR-risk_allele and adjust effect size sign accordingly 
oleate_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
oleate_kora_2$Effect2 <- ifelse(oleate_kora_2$Allele1 == oleate_kora_2$IR_allele,oleate_kora_2$Effect,-oleate_kora_2$Effect)
tyrosine_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
tyrosine_kora_2$Effect2 <- ifelse(tyrosine_kora_2$Allele1 == tyrosine_kora_2$IR_allele,tyrosine_kora_2$Effect,-tyrosine_kora_2$Effect)
palmit_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
palmit_kora_2$Effect2 <- ifelse(palmit_kora_2$Allele1 == palmit_kora_2$IR_allele,palmit_kora_2$Effect,-palmit_kora_2$Effect)
hippuric_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
hippuric_kora_2$Effect2 <- ifelse(hippuric_kora_2$Allele1 == hippuric_kora_2$IR_allele,hippuric_kora_2$Effect,-hippuric_kora_2$Effect)
mag182_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
mag182_kora_2$Effect2 <- ifelse(mag182_kora_2$Allele1 == mag182_kora_2$IR_allele,mag182_kora_2$Effect,-mag182_kora_2$Effect)
mag181_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
mag181_kora_2$Effect2 <- ifelse(mag181_kora_2$Allele1 == mag181_kora_2$IR_allele,mag181_kora_2$Effect,-mag181_kora_2$Effect)
gtoco_kora_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
gtoco_kora_2$Effect2 <- ifelse(gtoco_kora_2$Allele1 == gtoco_kora_2$IR_allele,gtoco_kora_2$Effect,-gtoco_kora_2$Effect)

# write.table(oleate_kora_2,"oleate_kora_2.txt"); # write.table(tyrosine_kora_2,"tyrosine_kora_2.txt")
# write.table(hippuric_kora_2,"hippuric_kora_2.txt"); # write.table(palmit_kora_2,"palmit_kora_2.txt")
# write.table(mag182_kora_2,"mag182_kora_2.txt"); # write.table(mag181_kora_2,"mag181_kora_2.txt")
# write.table(gtoco_kora_2,"gtoco_kora_2.txt")

  ## MR analysis with GRS based on summary-GWAS data with gtx package ####
  ## standardize to SD-unit
b_o_kora <- oleate_kora_2$Effect2 / (k_oleate_sd )
s_o_kora <- oleate_kora_2$StdErr / (k_oleate_sd )
b_t_kora <- tyrosine_kora_2$Effect2 / (k_tyrosine_sd )
s_t_kora <- tyrosine_kora_2$StdErr / (k_tyrosine_sd )
b_p_kora <- palmit_kora_2$Effect2 / (k_palmit_sd )
s_p_kora <- palmit_kora_2$StdErr / (k_palmit_sd )
b_h_kora <- hippuric_kora_2$Effect2 / (k_hippuric_sd )
s_h_kora <- hippuric_kora_2$StdErr / (k_hippuric_sd )
b_m182_kora <- mag182_kora_2$Effect2 / (k_mag182_sd )
s_m182_kora <- mag182_kora_2$StdErr / (k_mag182_sd )
b_m181_kora <- mag181_kora_2$Effect2 / (k_mag181_sd )
s_m181_kora <- mag181_kora_2$StdErr / (k_mag181_sd )
b_g_kora <- gtoco_kora_2$Effect2 / (k_gtoco_sd )
s_g_kora <- gtoco_kora_2$StdErr / (k_gtoco_sd )

n <- rep(1,10); wi <- rep(1,10)
o_kora <- data.frame(cbind(oleate_kora_2$MarkerName, b_o_kora,s_o_kora,n,wi))
t_kora <- data.frame(cbind(tyrosine_kora_2$MarkerName,b_t_kora,s_t_kora,n,wi))
p_kora <- data.frame(cbind(palmit_kora_2$MarkerName,b_p_kora,s_p_kora,n,wi))
h_kora <- data.frame(cbind(hippuric_kora_2$MarkerName,b_h_kora,s_h_kora,n,wi))
m182_kora <- data.frame(cbind(mag182_kora_2$MarkerName, b_m182_kora,s_m182_kora,n,wi))
m181_kora <- data.frame(cbind(mag181_kora_2$MarkerName, b_m181_kora,s_m181_kora,n,wi))
g_kora <- data.frame(cbind(gtoco_kora_2$MarkerName, b_g_kora,s_g_kora,n,wi))

attach(o_kora)
yko = grs.summary(wi, b_o_kora, s_o_kora, n)
grs.plot(wi, b_o_kora, s_o_kora, oleate_kora_2$MarkerName)
IR_MR_beta <- yko$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((yko$aSE/yko$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_oleate_kora <- cbind("kora oleate",yko$ahat, yko$aSE,yko$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(t_kora)
yto = grs.summary(wi, b_t_kora, s_t_kora, n)
grs.plot(wi, b_t_kora, s_t_kora, tyrosine_kora_2$MarkerName)
IR_MR_beta <- yto$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((yto$aSE/yto$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_tyrosine_kora <- cbind("kora tyrosine",yto$ahat, yto$aSE,yto$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(p_kora)
ypo = grs.summary(wi, b_p_kora, s_p_kora, n)
grs.plot(wi, b_p_kora, s_p_kora, palmit_kora_2$MarkerName)
IR_MR_beta <- ypo$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ypo$aSE/ypo$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_palmit_kora <- cbind("kora palmitoleate",ypo$ahat, ypo$aSE,ypo$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(h_kora)
yho = grs.summary(wi, b_h_kora, s_h_kora, n)
grs.plot(wi, b_h_kora, s_h_kora, hippuric_kora_2$MarkerName)
IR_MR_beta <- yho$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((yho$aSE/yho$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_hippuric_kora <- cbind("kora hippurate",yho$ahat, yho$aSE,yho$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(m182_kora)
ym182o = grs.summary(wi, b_m182_kora, s_m182_kora, n)
grs.plot(wi, b_m182_kora, s_m182_kora, mag182_kora_2$MarkerName)
IR_MR_beta <- ym182o$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ym182o$aSE/ym182o$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_m182_kora <- cbind("kora 1-monolinolein (close to MAG18_2)",ym182o$ahat, ym182o$aSE,ym182o$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(m181_kora)
ym181o = grs.summary(wi, b_m181_kora, s_m181_kora, n)
grs.plot(wi, b_m181_kora, s_m181_kora, mag181_kora_2$MarkerName)
IR_MR_beta <- ym181o$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ym181o$aSE/ym181o$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_m181_kora <- cbind("kora MAG18_1",ym181o$ahat, ym181o$aSE,ym181o$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(g_kora)
ygo = grs.summary(wi, b_g_kora, s_g_kora, n)
grs.plot(wi, b_g_kora, s_g_kora, gtoco_kora_2$MarkerName)
IR_MR_beta <- ygo$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ygo$aSE/ygo$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_gtoco_kora <- cbind("kora gamma-tocopherol",ygo$ahat, ygo$aSE,ygo$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

MR_kora <- rbind(MR_oleate_kora,MR_palmit_kora,MR_hippuric_kora,MR_tyrosine_kora,MR_m181_kora,MR_m182_kora,MR_gtoco_kora)
colnames(MR_kora) <- c("metabolte","GRS_beta","GRS_se","GRS_pval","IV_beta","IV_se","IR_pval")

# write.table(MR_kora,"KORA_IR_MR_results.txt")

#### CHARGE consortium, Wu et al., 161n7 = palmitoleate, 181n9 = oleate

palmit_charge <- read.delim("CHARGE_161n7.txt")
oleate_charge <- read.delim("CHARGE_181n9.txt")

oleate_charge_2 <- oleate_charge[which(as.character(oleate_charge$MarkerName) %in% as.character(ir_snps$rs_id)),]
oleate_charge_2 <- oleate_charge_2[order(as.character(oleate_charge_2$MarkerName)),] 
palmit_charge_2 <- palmit_charge[which(as.character(palmit_charge$MarkerName) %in% as.character(ir_snps$rs_id)),]
palmit_charge_2 <- palmit_charge_2[order(as.character(palmit_charge_2$MarkerName)),] 

oleate_charge_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
oleate_charge_2$Effect2 <- ifelse(oleate_charge_2$Allele1 == oleate_charge_2$IR_allele,oleate_charge_2$Effect,-oleate_charge_2$Effect)
palmit_charge_2$IR_allele <- (tolower(as.character(ir_snps$IR_allele)))
palmit_charge_2$Effect2 <- ifelse(palmit_charge_2$Allele1 == palmit_charge_2$IR_allele,palmit_charge_2$Effect,-palmit_charge_2$Effect)

# write.table(oleate_charge_2,"oleate_charge_2.txt")
# write.table(palmit_charge_2,"palmit_charge_2.txt")

  ## SD from Table 1 of Wu et al. 
c_oleate_sd <- sqrt(((1.1^2)*3269 + (1.2^2)*1507 + (1.1^2)*2404 + (3.7^2)*1075 + (1.2^2)*706) / (3269+1507+2404+1075+706))
c_palmit_sd <- sqrt(((0.18^2)*3269 + (0.22^2)*1507 + (0.21^2)*2404 + (0.9^2)*1075 + (0.25^2)*706) / (3269+1507+2404+1075+706))

b_o_charge <- oleate_charge_2$Effect2 / (c_oleate_sd )
s_o_charge <- oleate_charge_2$StdErr / (c_oleate_sd )
b_p_charge <- palmit_charge_2$Effect2 / (c_palmit_sd )
s_p_charge <- palmit_charge_2$StdErr / (c_palmit_sd )

o_charge <- data.frame(cbind(oleate_charge_2$MarkerName, b_o_charge,s_o_charge,n,wi))
p_charge <- data.frame(cbind(palmit_charge_2$MarkerName, b_p_charge,s_p_charge,n,wi))

attach(o_charge)
yco = grs.summary(wi, b_o_charge, s_o_charge, n)
grs.plot(wi, b_o_charge, s_o_charge, oleate_charge_2$MarkerName)
IR_MR_beta <- yco$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((yco$aSE/yco$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_oleate_charge <- cbind("charge oleate",yco$ahat, yco$aSE,yco$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

attach(p_charge)
ycp = grs.summary(wi, b_p_charge, s_p_charge, n)
grs.plot(wi, b_p_charge, s_p_charge, palmit_charge_2$MarkerName)
IR_MR_beta <- ycp$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ycp$aSE/ycp$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_palmit_charge <- cbind("charge palmit",ycp$ahat, ycp$aSE,ycp$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

MR_charge <- rbind(MR_palmit_charge,MR_oleate_charge)
colnames(MR_charge) <- c("metabolte","GRS_beta","GRS_se","GRS_pval","IV_beta","IV_se","IR_pval")

# write.table(MR_charge,"IR_CHARGE_results.txt")

####  Finnish Tyrosine - GWAS results extracted by Ripatti team - adj. for age + sex; z-standardized 
  ## CAVE: effect allele is allele_B

tyrosine_fin <- read.delim("tukiainen_ir_tyrosine.txt")
tyrosine_fin_2 <- tyrosine_fin[order(as.character(tyrosine_fin$rsid)),] 
tyrosine_fin_2$IR_allele <- (as.character(ir_snps$IR_allele))
tyrosine_fin_2$Effect2 <- ifelse(tyrosine_fin_2$allele_B == tyrosine_fin_2$IR_allele,tyrosine_fin_2$beta,-tyrosine_fin_2$beta)

# write.table(tyrosine_fin_2,"tyrosine_fin_2.txt")

b_t_fin <- tyrosine_fin_2$Effect2
s_t_fin <- tyrosine_fin_2$se
t_fin <- data.frame(cbind(tyrosine_fin_2$rsid, b_t_fin,s_t_fin,n,wi))

attach(t_fin)
ytf = grs.summary(wi, b_t_fin, s_t_fin, n)
grs.plot(wi, b_t_fin, s_t_fin, tyrosine_fin_2$rsid)
IR_MR_beta <- ytf$ahat / SCOTTclamp_beta
IR_MR_se <- abs(IR_MR_beta)*sqrt( ((ytf$aSE/ytf$ahat)^2)+ ((SCOTTclamp_se/SCOTTclamp_beta)^2) )
IR_MR_z <- IR_MR_beta / IR_MR_se
IR_MR_p <- 2*pnorm(-abs(IR_MR_z))
MR_tyrosine_fin <- cbind("fin tyrosine",ytf$ahat, ytf$aSE,ytf$pval,IR_MR_beta,IR_MR_se,IR_MR_p)

# write.table(MR_tyrosine_fin,"IR_FINNISH_results.txt")

##### FOREST PLOTS 

### order
## palmit     SWE 2   KORA 2  CHARGE 1
## oleic      SWE 5   KORA 1  CHARGE 2
## hippuric   SWE 1   KORA 3
## tyrosine   SWE 3   KORA 4  FINNish
## MAG181     SWE 8   KORA 5
## MAG182     SWE 4   KORA 6
## gammaToco  SWE 9   KORA 7
## MAG140     SWE 6
## 3a0 Trep   SWE 7


temp3 <- read.table("IR_MR_U70_P70_TWGE_results.txt")
MR_charge <-  read.table("IR_CHARGE_results.txt")
MR_tyrosine_fin <- read.table("IR_FINNISH_results.txt")
MR_kora <- read.table("KORA_IR_MR_results.txt")
temp3 <- data.frame(temp3)
SWE <- subset(temp3,as.numeric(as.character(temp3$IR_MR_p)) < 0.1)
SWE$IR_MR_beta <- as.numeric(as.character(SWE$IR_MR_beta))
SWE$IR_MR_se <- as.numeric(as.character(SWE$IR_MR_se))
KORA <- data.frame(MR_kora)
KORA$IV_beta <- as.numeric(as.character(KORA$IV_beta))
KORA$IV_se <- as.numeric(as.character(KORA$IV_se))
CHARGE <- data.frame(MR_charge)
CHARGE$IV_beta <- as.numeric(as.character(CHARGE$IV_beta))
CHARGE$IV_se <- as.numeric(as.character(CHARGE$IV_se))
  
  ## upload insulin secretion results for plot
isSWE<- data.frame(read.table("Insecr_MR_results.txt"))
isKORA <- data.frame(read.table("KORA_IS_MR_results.txt"))

x <-c(SWE$IR_MR_beta[2],KORA$IV_beta[2],CHARGE$IV_beta[1],
      SWE$IR_MR_beta[5],KORA$IV_beta[1],CHARGE$IV_beta[2],
      SWE$IR_MR_beta[1],KORA$IV_beta[3],
      SWE$IR_MR_beta[3],KORA$IV_beta[4],MR_tyrosine_fin[,5],
      SWE$IR_MR_beta[8],KORA$IV_beta[5],
      SWE$IR_MR_beta[4],KORA$IV_beta[6],
      SWE$IR_MR_beta[9],KORA$IV_beta[7],
      SWE$IR_MR_beta[6],
      SWE$IR_MR_beta[7],
      isSWE[which(isSWE$X=="Trihydroxy_5b_cholanoic_acid"),]$Insecr_MR_beta, 
      isSWE[which(isSWE$X=="Bilirubin"),]$Insecr_MR_beta, isKORA$IV_beta[3])
sei <-c(SWE$IR_MR_se[2],KORA$IV_se[2],CHARGE$IV_se[1],
        SWE$IR_MR_se[5],KORA$IV_se[1],CHARGE$IV_se[1],
        SWE$IR_MR_se[1],KORA$IV_se[3],
        SWE$IR_MR_se[3],KORA$IV_se[4],MR_tyrosine_fin[,6],
        SWE$IR_MR_se[8],KORA$IV_se[5],
        SWE$IR_MR_se[4],KORA$IV_se[6],
        SWE$IR_MR_se[9],KORA$IV_se[7],
        SWE$IR_MR_se[6],
        SWE$IR_MR_se[7],
        isSWE[which(isSWE$X=="Trihydroxy_5b_cholanoic_acid"),]$Insecr_MR_se,
        isSWE[which(isSWE$X=="Bilirubin"),]$Insecr_MR_se, isKORA$IV_se[3])
M <- c("Palmitoleic acid","","",
       "Oleic acid","","",
       "Hippuric acid","",
       "Tyrosine","","",
       "MG(18:1)","",
       "MG(18:2)","",
       "Gamma-tocopherol","",
       "MG(14:0)","Trihydroxy-cholanoic acid",
       "Trihydroxy-cholanoic acid",
       "Bilirubin","")
Study <- c(rep(c("Swedish Cohorts","KORA/TwinsUK","CHARGE"),2),c("Swedish Cohorts","KORA/TwinsUK"),
           c("Swedish Cohorts","KORA/TwinsUK","Finnish Cohorts"),rep(c("Swedish Cohorts","KORA/TwinsUK"),3),
           rep("Swedish Cohorts",3),c("Swedish Cohorts", "KORA/TwinsUK"))

# 10.71 8.95 pdf 13.5x12
par(mfrow=c(1,1))
forest.default(x,sei=sei,psize=0.85,slab=M,
               ilab=Study,ilab.xpos=-2.6,cex=0.85,
               xlab="Causal estimate (95% CI) in SD-units",
               ylim=c(0,46.5),
               rows=c(41,40,39,  36,35,34,   31,30,   27,26,25,   22,21,  18,17,  14,13,
                      10,  8,  4,  2,1))
addpoly.default(rma.uni(yi = x[1:3],sei = sei[1:3])$b, sei=rma.uni(yi = x[1:3],sei = sei[1:3])$se, mlab=" ",rows=38,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[4:6],sei = sei[4:6])$b, sei=rma.uni(yi = x[4:6],sei = sei[4:6])$se, mlab=" ",rows=33,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[7:8],sei = sei[7:8])$b, sei=rma.uni(yi = x[7:8],sei = sei[7:8])$se, mlab=" ",rows=29,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[9:11],sei = sei[9:11])$b, sei=rma.uni(yi = x[9:11],sei = sei[9:11])$se, mlab=" ",rows=24,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[12:13],sei = sei[12:13])$b, sei=rma.uni(yi = x[12:13],sei = sei[12:13])$se, mlab=" ",rows=20,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[14:15],sei = sei[14:15])$b, sei=rma.uni(yi = x[14:15],sei = sei[14:15])$se, mlab=" ",rows=16,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[16:17],sei = sei[16:17])$b, sei=rma.uni(yi = x[16:17],sei = sei[16:17])$se, mlab=" ",rows=12,cex=0.85,font=2,col="grey")
addpoly.default(rma.uni(yi = x[21:22],sei = sei[21:22])$b, sei=rma.uni(yi = x[21:22],sei = sei[21:22])$se, mlab=" ",rows=0,cex=0.85,font=2,col="grey")
text(-4.4,43,"Insulin Resistance",cex=1,font=2)
text(-4.09,6,"Impaired Insulin Secretion",cex=1,font=2)


