rm(list=ls())

# Rfit analyses

library(dplyr)
library(ggplot2)

#CAG.RA <- read.table("bali_abundance_byCAG.csv",row.names = 1, header = T,sep = ",")
CAG.RA <- read.table("data/bali_coabundance_RA_byCAG.csv",row.names = 1, header = T,sep = ",")

head(CAG.RA)

dietCPF <- read.table("data/bali_food_recall_all_summary(n=41).csv",sep = ",",header = T)
head(dietCPF)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

dietCPF %>%
  mutate(outlier = ifelse(is_outlier(energy.kcal), as.character(SampleID), NA)) %>%
  ggplot(., aes(x = "Energy Intake (kcal)", y = energy.kcal)) +
  theme_bw()+
  geom_boxplot(width=0.8) +
  geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) +
  xlab("") + ylab("")

apply(dietCPF[,c("energy.kcal","carbs.percentage","fat.percentage","prot.percentage")],2, function(x) shapiro.test(x)) -> shapiout
do.call(rbind,lapply(shapiout, function(x) unlist(x[c("statistic","p.value")]))) -> shapiout
shapiout

library(psych)
corr.test(dietCPF[,c("energy.kcal","carbs.percentage","fat.percentage","prot.percentage")],method = "kendall",adjust = "fdr")

dietCPF$SampleID <- as.character(dietCPF$SampleID)
n <- dietCPF$SampleID[is_outlier(dietCPF$energy.kcal)]
dietCPF <- dietCPF[-which(dietCPF$SampleID == n),]

# Load cluster data
unifrac.clust <- read.csv("bali_cluster_ward's_rare90k.csv")
unifrac.clust$Order <- NULL
unifrac.clust$Cluster <- factor(unifrac.clust$Cluster, c("Type-B","Type-P"))

rfit.combine <- merge.data.frame(unifrac.clust, dietCPF[,c("SampleID","energy.kcal","carbs.percentage","fat.percentage","prot.percentage")], by="SampleID")

# Load BMI Data
antrop <- read.csv("data/bali_antropometry.csv",header = T,sep = ",")
antrop <- antrop[-which(antrop$SampleID == n),c("SampleID","BMI")]

rfit.combine <- merge.data.frame(rfit.combine, antrop, by="SampleID")

head(rfit.combine)
str(rfit.combine)


library(Rfit)
library(ggplot2)
library(gtools)

# calculate adjusted critical values for 95% confidence intervals
cv2 <- confintadjust(40,2,alpha=0.05, method="bonferroni")$cv #where sample size (n) = 40, and number of comparisons (k) = 5 for carb, fat, prot, bmi, and comm.type
cv3 <- confintadjust(40,3,alpha=0.05, method="bonferroni")$cv #where sample size (n) = 40, and number of comparisons (k) = 5 for carb, fat, prot, bmi, and comm.type
cv4 <- confintadjust(40,4,alpha=0.05, method="bonferroni")$cv #where sample size (n) = 40, and number of comparisons (k) = 5 for carb, fat, prot, bmi, and comm.type
cv5 <- confintadjust(40,5,alpha=0.05, method="bonferroni")$cv #where sample size (n) = 40, and number of comparisons (k) = 5 for carb, fat, prot, bmi, and comm.type
cv6 <- confintadjust(40,6,alpha=0.05, method="bonferroni")$cv #where sample size (n) = 40, and number of comparisons (k) = 5 for carb, fat, prot, bmi, and comm.type

cv <- c(cv2, cv3, cv4, cv5, cv6)
names(cv) <- c("k2","k3","k4","k5","k6")
print(cv)

fit.coef.all <- data.frame(response=c(), 
                           predictors=c(), 
                           Estimate=c(), 
                           Std..Error=c(), 
                           t.value=c(), 
                           p.value=c(), 
                           signif=c(), 
                           LC=c(), 
                           UC=c(), 
                           cols=c(), 
                           model=c(),
                           model.R2.multi=c(), 
                           model.drop.pval=c())

for(i in colnames(CAG.RA)){
  rfit.data <- CAG.RA[i]
  colnames(rfit.data) <- "response"
  rfit.data$SampleID <- rownames(rfit.data)
  rfit.data <- merge.data.frame(rfit.data, rfit.combine, by="SampleID")
  head(rfit.data)
  
  rownames(rfit.data) <- rfit.data$SampleID
  rfit.data$SampleID <- NULL
  
  print(paste("response =",i))
  
  fit1 <- rfit(response ~  carbs.percentage + fat.percentage + prot.percentage, data = rfit.data) 
  fit1.res <- summary(fit1,overall.test = "drop")
  print(fit1.res)
  
  fit2 <- rfit(response ~ BMI + carbs.percentage  + fat.percentage + prot.percentage, data = rfit.data) 
  fit2.res <- summary(fit2,overall.test = "drop")
  print(fit2.res)
  
  fit3 <- rfit(response ~ Cluster + carbs.percentage  + fat.percentage + prot.percentage, data = rfit.data) 
  fit3.res <- summary(fit3,overall.test = "drop")
  print(fit3.res)
  
  fit4 <- rfit(response ~ Cluster + BMI + carbs.percentage + fat.percentage + prot.percentage, data = rfit.data) 
  fit4.res <- summary(fit4,overall.test = "drop")
  print(fit4.res)
  
  fit.ls <- list(fit1 = data.frame(fit1.res$coefficients, model.R2.multi=fit1.res$R2, model.drop.pval=fit1.res$droppval), 
                 fit2 = data.frame(fit2.res$coefficients, model.R2.multi=fit2.res$R2, model.drop.pval=fit2.res$droppval),
                 fit3 = data.frame(fit3.res$coefficients, model.R2.multi=fit3.res$R2, model.drop.pval=fit3.res$droppval), 
                 fit4 = data.frame(fit4.res$coefficients, model.R2.multi=fit4.res$R2, model.drop.pval=fit4.res$droppval))
  
  lab <- paste(i, "~",c("Diet", "Diet + BMI", "Diet + Comm.Type", "Diet + BMI + Comm.Type"))
  names(fit.ls) <- lab

  for(n in c(1:4)){
      fit.coef <- fit.ls[[names(fit.ls)[n]]]
      fit.coef <- data.frame(response=i, predictors=rownames(fit.coef), fit.coef)
      rownames(fit.coef) <- NULL
      fit.coef$predictors <- gsub("Cluster","",fit.coef$predictors)
      fit.coef$predictors <- gsub(".percentage", "%",fit.coef$predictors)
      fit.coef$predictors <- factor(fit.coef$predictors, c("(Intercept)","Type-P","BMI","carbs%","fat%","prot%"))
      
      fit.coef$signif <- as.character(stars.pval(fit.coef$p.value))
      
      k = paste("k",nrow(fit.coef)-1, sep="")
      fit.coef$LC <- fit.coef$Estimate-(cv[k]*fit.coef$Std..Error)
      fit.coef$UC <- fit.coef$Estimate+(cv[k]*fit.coef$Std..Error)

      fit.coef$cols <- ifelse(fit.coef$Estimate < 0, "red","royalblue")
      
      fit.coef$model <- names(fit.ls)[n]
      
      fit.coef <- fit.coef[,c(1:6,9,10:13,7:8)]
      
      fit.coef.all <- rbind(fit.coef.all, fit.coef)
    }
}

print(fit.coef.all)
tail(fit.coef.all)

#write.table(fit.coef.all, "Figures/bali_CAG_Rfit_Collate_Result_Macronutrient+BMI+Comm.Type.csv",quote = F,sep = ",",row.names = F)

# Plot result
fit.coef.all$facet.col <- unlist(lapply(strsplit(as.character(fit.coef.all$model)," ~ "), function(x) x[2]))
fit.coef.all$facet.col <- factor(fit.coef.all$facet.col)
fit.coef.all$facet.col <- factor(fit.coef.all$facet.col, levels(fit.coef.all$facet.col) [c(1,4,2,3)])

fit.coef.all$model <- factor(fit.coef.all$model)
fit.coef.all$response <- factor(fit.coef.all$response)

fit.coef.all$label <- paste(as.character(round(fit.coef.all$Estimate,2)),
                            gsub(pattern = "\\.",replacement = "m",x = as.character(stars.pval(fit.coef.all$p.value))),
                        sep = "")


fit.coef.all$model.R2.multi <- round(fit.coef.all$model.R2.multi, 2)
fit.coef.all$model.drop.pval <- round(fit.coef.all$model.drop.pval, 3)

fit.coef.all$predictors <- factor(fit.coef.all$predictors, c("(Intercept)","Type-P","BMI","carbs%","fat%","prot%"))

it <- which(fit.coef.all$predictors == "(Intercept)")

fit.subs <- fit.coef.all[-it,]
#fit.subs <- fit.subs[fit.subs$response %in% c("CAG1","CAG2","CAG3","CAG4","CAG7"),]

ggplot(fit.subs, aes(y=Estimate, x=predictors,col=cols)) +
  theme_bw() +
  facet_grid(response~facet.col) +
  geom_hline(yintercept = 0, col="grey85", size=2) +
  scale_x_discrete(limits = rev(levels(fit.coef.all$predictors)[-1]))+
  scale_y_continuous(limits = c(-30,40),breaks=sort(c(0,seq(-40,40,by = 20))))+
  geom_errorbar(aes(ymin=LC, ymax=UC), width=0, size=1.2) + 
  geom_point(size=3, shape=1) + 
  scale_color_manual(values = c(red="red",royalblue="royalblue"))+
  geom_text(aes(label=label), col="black",nudge_x = 0.3, size=3) +
  coord_flip() +
  geom_text(aes(label=paste0("R2=",model.R2.multi,"\n p=",model.drop.pval)),x = 2, y=30, col="black", size=3)+
  theme(legend.position = "none",
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank())

