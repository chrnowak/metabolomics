##########################################################################################
#### Impaired insulin secretion ##########################################################

dir <- (" ... /Scripts")
setwd(dir)

library(gtx)

metabonames <- c("MAG_18_2","MAG_16_1","MAG_18_1","Piperine","Caffeine","LPE_18_2","DOCA_Glycine_Conj","LPE_18_1","MAG_16_0",
                 "Corticosterone","LPC_18_3","LPC_20_0","gamma_Glu_Leu","MAG_14_0","LPC_18_2_0_0","Propanolol","LPC_0_0_18_2",
                 "DAG_34_1","Arachidonic_acid_ethyl_ester","PC_36_4a","Eicosatrienoic_acid","LPE_18_1b","PC_16_2","gamma_Tocopherol",
                 "Cinnamic_acid","LPC_18_1a","LPC_18_1b","Paraxanthine_Theophylline","L_Tyrosine","Hippuric_acid","linoleoyl_2_stearoyl_sn_glycerol",
                 "Chenodeoxycholic_acid_Glyc_Conj","Decanoyl_L_carnitine","LPE_20_4","LPE_18_1_P","Glycocholic_acid","L_Acetylcarnitine","Betaine",
                 "Trihydroxy_5b_cholanoic_acid","PC_38_7","PC_38_3","Oleic_acid","LPE_16_0a","Choline","LPC_20_3a","Treprostinil","Hexose","Docosapentaenoic_acid",
                 "Bilirubin","Bilirubin_IXa","Palmitoleic_acid","Limoleic_acid")

COMB <- read.table("P_TWGE_nonDM_for_MR.txt")
varnames <- c(metabonames,"SEX","AGE","Cohort_1","CNT2")
COMB2 <- COMB[,which(names(COMB) %in% varnames)]
is_snps <- read.delim("scott_is_snps.txt")

  ## Scott et al. GRS score association with insulinogenic index 
  ## sign reversed here, so that increased GRS = impaired insulin secretion
SCOTTinsecr_beta <- 0.05 
SCOTTinsecr_se <- (-0.06 - -0.04)/3.92

### requested by reviewers: population stratification adjustment
### adjustment for first 3 genetic principle components

pivus_PC <- read.table("~/Desktop/pivus.+.qced.pca.evec.txt", quote="\"")
twge_PC <- data.frame(read.table("~/Desktop/twge.+.pca.evec.txt", quote="\""))
split_1 <- str_split_fixed(as.character(twge_PC[,1]), ":", 2)
twge_PC$FID <- (split_1[,1])
pivus_PC$lpnr <- (as.numeric(pivus_PC[,1]))

COMB2$PC1<-c()
COMB2$PC2<-c()
COMB2$PC3<-c()
for (i in 1:nrow(COMB2)){
  if (!is.na(COMB2$lpnr[i])) COMB2$PC1[i] <- pivus_PC[match(COMB2$lpnr[i],pivus_PC$lpnr),2]
  else COMB2$PC1[i] <- twge_PC[match(COMB2$gwas_fid[i],twge_PC$FID),2]
  if (!is.na(COMB2$lpnr[i])) COMB2$PC2[i] <- pivus_PC[match(COMB2$lpnr[i],pivus_PC$lpnr),3]
  else COMB2$PC2[i] <- twge_PC[match(COMB2$gwas_fid[i],twge_PC$FID),3]
  if (!is.na(COMB2$lpnr[i])) COMB2$PC3[i] <- pivus_PC[match(COMB2$lpnr[i],pivus_PC$lpnr),4]
  else COMB2$PC3[i] <- twge_PC[match(COMB2$gwas_fid[i],twge_PC$FID),4]
}

COMB2$PC1 <- ifelse(is.na(COMB2$PC1),mean(COMB2$PC1,na.rm=T),COMB2$PC1)
COMB2$PC2 <- ifelse(is.na(COMB2$PC2),mean(COMB2$PC2,na.rm=T),COMB2$PC2)
COMB2$PC3 <- ifelse(is.na(COMB2$PC3),mean(COMB2$PC3,na.rm=T),COMB2$PC3)

CNT2_Insecr_Beta <- sapply(which(names(COMB2) %in% metabonames), function(x) 
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1  + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$coef[5,1])
CNT2_Insecr_SE <- sapply(which(names(COMB2) %in% metabonames), function(x) 
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1  + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$coef[5,2])
CNT2_Insecr_P <- sapply(which(names(COMB2) %in% metabonames), function(x) 
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1  + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$coef[5,4])
CNT2_Insecr_Fstat <- sapply(which(names(COMB2) %in% metabonames), function(x) 
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1 + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$fstatistic[1])
CNT2_Insecr_Fstat_P <- sapply(which(names(COMB2) %in% metabonames), function(x) pf(
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1 + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$fstatistic[1],
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1 + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$fstatistic[2],
  summary(lm(scale(COMB2[,x]) ~ COMB2$AGE + COMB2$SEX + COMB2$Cohort_1 + COMB2$CNT2 + COMB2$PC1 + COMB2$PC2 + COMB2$PC3))$fstatistic[3],
  lower.tail=F)) # model F-statistic extraction
CNT2_Insec_results <- cbind(names(COMB2)[which(names(COMB2) %in% metabonames)], CNT2_Insecr_Beta,CNT2_Insecr_SE,CNT2_Insecr_P,CNT2_Insecr_Fstat,CNT2_Insecr_Fstat_P)

# write.table(CNT2_Insec_results,"IS_metabo_GRS_assoc.txt",quote=F)

Insecr_MR_beta <- sapply(1:length(metabonames), function(x) CNT2_Insecr_Beta[x] / SCOTTinsecr_beta)
Insecr_MR_se <- sapply(1:length(metabonames), function(x) 
  abs(Insecr_MR_beta[x])*sqrt( ((CNT2_Insecr_SE[x]/CNT2_Insecr_Beta[x])^2)+ ((SCOTTinsecr_se/SCOTTinsecr_beta)^2) ))
Insecr_MR_z <- sapply(1:length(metabonames), function(x) Insecr_MR_beta[x] / Insecr_MR_se[x])
Insecr_MR_p <- sapply(1:length(metabonames), function(x) 2*pnorm(-abs(Insecr_MR_z[x])))
Insecr_MR_results <- cbind(names(COMB2)[which(names(COMB2) %in% metabonames)],Insecr_MR_beta,Insecr_MR_se,Insecr_MR_z,Insecr_MR_p)
Insecr_MR_results <- Insecr_MR_results[order(Insecr_MR_results[,5]),]

# write.table(Insecr_MR_results, "Insecr_MR_results.txt")

##### KORA/TwinsUK: Bilirubin (Z,Z): M27716; SDs: 0.320 (n=5092)  0.312 (n=1720)

B3_kora <- read.delim("M27716.metal.pos.txt")
k_b3_sd <- sqrt(((0.320^2)*5092 + (0.312^2)*1720)/(5092+1720))
is_snps <- data.frame(read.delim("scott_is_snps.txt")) # Scott et al. Suppl Table
is_snps <- is_snps[order(as.character(is_snps$rs_id)),] 
B3_kora_2 <- B3_kora[which(as.character(B3_kora$MarkerName) %in% as.character(is_snps$rs_id)),]
B3_kora_2 <- B3_kora_2[order(as.character(B3_kora_2$MarkerName)),] # order alph by rs_id

# write.table(B3_kora_2,"B3_kora_2.txt")

  ## snp ids in KORA
is_snps2 <- is_snps[which(as.character(is_snps$rs_id) %in% as.character(B3_kora$MarkerName)),]
B3_kora_2$IR_allele <- (tolower(as.character(is_snps2$IR_allele)))
B3_kora_2$Effect2 <- ifelse(B3_kora_2$Allele1 == B3_kora_2$IR_allele,B3_kora_2$Effect,-B3_kora_2$Effect)

n <- rep(1,nrow(is_snps2)); wi <- rep(1,nrow(is_snps2))
b_b3_kora <- B3_kora_2$Effect2 / (k_b3_sd )
s_b3_kora <- B3_kora_2$StdErr / (k_b3_sd )
b3_kora <- data.frame(cbind(as.character(B3_kora_2$MarkerName), b_b3_kora,s_b3_kora,n,wi))

# write.table(b3_kora,"b3_kora.txt")

attach(b3_kora)
y3 = grs.summary(wi, b_b3_kora, s_b3_kora, n)
grs.plot(wi, b_b3_kora, s_b3_kora, B3_kora_2$MarkerName)
is_MR_beta <- y3$ahat / SCOTTinsecr_beta
is_MR_se <- abs(is_MR_beta)*sqrt( ((y3$aSE/y3$ahat)^2)+ ((SCOTTinsecr_se/SCOTTinsecr_beta)^2) )
is_MR_z <- is_MR_beta / is_MR_se
is_MR_p <- 2*pnorm(-abs(is_MR_z))
MR_B3_kora <- cbind("kora bilirubin Z,Z",y3$ahat, y3$aSE,y3$pval,is_MR_beta,is_MR_se,is_MR_p)

tab0 <- (MR_B3_kora)
colnames(tab0) <- c("metabolite","GRS_beta","GRS_se","GRS_pval","IV_beta","IV_se","IR_pval")

# write.table(tab0,"KORA_IS_MR_results.txt")

