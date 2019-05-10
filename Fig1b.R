#### Beta Diversity analysis of 6 populations (Fig. 1b)

library(dplyr)
library(ggplot2)
library(reshape2)
library(vegan)

#Load environmental data
environ.var <- read.csv("data/6_populations_metadata_tidy.csv", header=TRUE, row.names=1,stringsAsFactors = T)

# setting factorial variables
popord <- c("Balinese","Hadza","Malawian","Guahibo","Italian","American")
environ.var$Population <- factor(environ.var$Population, popord)
environ.var$Obesity <- ifelse(is.na(environ.var$BMI), "NA", ifelse(environ.var$BMI > 30, "obese","lean"))
environ.var$Obesity <- factor(environ.var$Obesity, c("lean","obese","NA"))
environ.var$Lifestyle <- factor(environ.var$Lifestyle, c("Transitional","Industrial","PreIndustrial"))
environ.var$Sex <- factor(environ.var$Sex, c("Male","Female"))
environ.var$femcount <- as.numeric(environ.var$Sex) -1

# set color by lifestyle (+ Bali)
lifestyle_colour <- c("red","lightskyblue","chartreuse3")
names(lifestyle_colour) <- levels(environ.var$Lifestyle)

# set color by country
country_colour <- c("red","deepskyblue","mediumspringgreen","gold","royalblue3","forestgreen")
names(country_colour) <- levels(environ.var$Country)

#Load distance matrix
w.unifrac.dm <- read.table("data/weighted_unifrac_rarefaction_4000_adults.txt")

reorder_dm <- function(var){
  var[order(rownames(var)),order(colnames(var))]
}

w.unifrac.dm <- reorder_dm(w.unifrac.dm)
w.unifrac.dist <- as.dist(w.unifrac.dm)

# match sample order in the dataset
sample_name <- rownames(w.unifrac.dm)
environ.var <- environ.var[sample_name,]

# create an ordination object
cap <- capscale(w.unifrac.dist ~ 1, data=environ.var)

# extract eigenvalues
eig <- cap$CA$eig
axisVar <- paste(round((eig/sum(eig)) * 100,1),"%")
names(axisVar) <- names(eig)

axis_lab <- c()
for(i in 1:(length(axisVar))){
  axis_lab <- c(axis_lab,paste(names(axisVar)[i]," (",axisVar[i],") variance explained", sep=""))
}

axis_lab <- gsub(pattern = "MDS",replacement = "Axis ",axis_lab)

# plot with color
par(mar=c(5,5,5,5))
plot(cap, type="n",display="sites", cex.axis=1.5, cex.lab=1.5, xlab=axis_lab[1], ylab=axis_lab[2])
for(i in 1:3){ordihull(cap, environ.var$Lifestyle, display = "sites", draw = "polygon", show.groups=levels(environ.var$Lifestyle)[i], col=lifestyle_colour[i],border = F, alpha=0.3)}
points(cap, disp="sites", pch=c(22,24,21)[environ.var$Lifestyle], bg=country_colour[environ.var$Country], cex=1.5)
legend("topright",pch=c(22,24,21,21,24,21),pt.bg = country_colour, legend = levels(environ.var$Country),xjust = 1,yjust=1, cex=1.2)
title("Fig. 1b - Weighted Unifrac PC1 vs PC2")


# Within Group Analysis
tmp <- data.frame(w.unifrac.dm,Sample1=rownames(w.unifrac.dm))
tmp2 <- data.frame(Sample1 = rownames(environ.var),Population1=environ.var$Population)

tmp <- merge.data.frame(tmp,tmp2, by="Sample1")
tmp <- melt(tmp, id.vars = c("Sample1","Population1"),variable.name = "Sample2",value.name = "w.UniFrac")
tmp2 <- data.frame(Sample2 = rownames(environ.var),Population2=environ.var$Population)
tmp <- merge.data.frame(tmp,tmp2, by="Sample2")
tmp <- tmp[,c(2,1,3,5,4)]
head(tmp)

withingroup <- tmp[tmp$Population1==tmp$Population2,]

p <- ggplot(withingroup, aes(x=Population1,y=w.UniFrac)) + theme_bw() +
  geom_violin(fill="grey80", col=F) + geom_boxplot(width=0.3,outlier.shape = 21)+
  xlab("Population") + ylab("Weighted UniFrac Distance")+ggtitle("Within-Population Distance")

print(p)

kruskal.test(x=withingroup$w.UniFrac, g = withingroup$Population1)
pairwise.wilcox.test(x=withingroup$w.UniFrac, g = withingroup$Population1,paired = F, p.adjust.method = "fdr") -> pairw
pairw <- data.frame(pairw$p.value)
pairw$Population1 <- rownames(pairw)
pairw <- melt(pairw, id.vars = "Population1",value.name = "pvalue",variable.name = "Population2")
pairw <- na.omit(pairw)

pairw$signif <- ifelse(pairw$pvalue <0.001, "***",
                       ifelse(pairw$pvalue <0.01, "**",
                              ifelse(pairw$pvalue <0.05, "*",
                                     "ns")))
print(pairw)

df.annot <- pairw
rownames(df.annot) <- NULL
df.annot$Population2 <- factor(df.annot$Population2, levels = c("Balinese","Hadza","Malawian","Guahibo","Italian","American"))
df.annot$Population1 <- factor(df.annot$Population1, levels = c("Balinese","Hadza","Malawian","Guahibo","Italian","American"))

df.annot <- df.annot[order(df.annot$Population2),]

increase <- 0.025
df.annot$x <- (as.numeric(df.annot$Population2) + as.numeric(df.annot$Population1))/2
df.annot$y <- 0.3+seq(from=0.05, to=(increase*(1+nrow(df.annot))), by=increase)

df.annot$ymin <- df.annot$y-0.005
df.annot$ymax <- df.annot$y-0.005

df.annot$xmin <- as.numeric(df.annot$Population2)
df.annot$xmax <- as.numeric(df.annot$Population1)

print(df.annot)

p + 
  annotate("segment",x = df.annot$xmin[1], xend = df.annot$xmax[1], y = df.annot$ymin[1], yend=df.annot$ymax[1])+
  annotate("segment",x = df.annot$xmin[2], xend = df.annot$xmax[2], y = df.annot$ymin[2], yend=df.annot$ymax[2])+
  annotate("segment",x = df.annot$xmin[3], xend = df.annot$xmax[3], y = df.annot$ymin[3], yend=df.annot$ymax[3])+
  annotate("segment",x = df.annot$xmin[4], xend = df.annot$xmax[4], y = df.annot$ymin[4], yend=df.annot$ymax[4])+
  annotate("segment",x = df.annot$xmin[5], xend = df.annot$xmax[5], y = df.annot$ymin[5], yend=df.annot$ymax[5])+
  annotate("segment",x = df.annot$xmin[6], xend = df.annot$xmax[6], y = df.annot$ymin[6], yend=df.annot$ymax[6])+
  annotate("segment",x = df.annot$xmin[7], xend = df.annot$xmax[7], y = df.annot$ymin[7], yend=df.annot$ymax[7])+
  annotate("segment",x = df.annot$xmin[8], xend = df.annot$xmax[8], y = df.annot$ymin[8], yend=df.annot$ymax[8])+
  annotate("segment",x = df.annot$xmin[9], xend = df.annot$xmax[9], y = df.annot$ymin[9], yend=df.annot$ymax[9])+
  annotate("segment",x = df.annot$xmin[10], xend = df.annot$xmax[10], y = df.annot$ymin[10], yend=df.annot$ymax[10])+
  annotate("segment",x = df.annot$xmin[11], xend = df.annot$xmax[11], y = df.annot$ymin[11], yend=df.annot$ymax[11])+
  annotate("segment",x = df.annot$xmin[12], xend = df.annot$xmax[12], y = df.annot$ymin[12], yend=df.annot$ymax[12])+
  annotate("segment",x = df.annot$xmin[13], xend = df.annot$xmax[13], y = df.annot$ymin[13], yend=df.annot$ymax[13])+
  annotate("segment",x = df.annot$xmin[14], xend = df.annot$xmax[14], y = df.annot$ymin[14], yend=df.annot$ymax[14])+
  annotate("segment",x = df.annot$xmin[15], xend = df.annot$xmax[15], y = df.annot$ymin[15], yend=df.annot$ymax[15])+
  annotate("segment",x=df.annot$xmin[1], xend=df.annot$xmin[1], y=df.annot$ymin[1]-0.01,yend = df.annot$ymin[1]) + annotate("segment", x=df.annot$xmax[1], xend=df.annot$xmax[1], y=df.annot$ymin[1]-0.01,yend = df.annot$ymin[1]) +
  annotate("segment",x=df.annot$xmin[2], xend=df.annot$xmin[2], y=df.annot$ymin[2]-0.01,yend = df.annot$ymin[2]) + annotate("segment", x=df.annot$xmax[2], xend=df.annot$xmax[2], y=df.annot$ymin[2]-0.01,yend = df.annot$ymin[2]) +
  annotate("segment",x=df.annot$xmin[3], xend=df.annot$xmin[3], y=df.annot$ymin[3]-0.01,yend = df.annot$ymin[3]) + annotate("segment", x=df.annot$xmax[3], xend=df.annot$xmax[3], y=df.annot$ymin[3]-0.01,yend = df.annot$ymin[3]) +
  annotate("segment",x=df.annot$xmin[4], xend=df.annot$xmin[4], y=df.annot$ymin[4]-0.01,yend = df.annot$ymin[4]) + annotate("segment", x=df.annot$xmax[4], xend=df.annot$xmax[4],y=df.annot$ymin[4]-0.01,yend = df.annot$ymin[4]) +
  annotate("segment",x=df.annot$xmin[5], xend=df.annot$xmin[5], y=df.annot$ymin[5]-0.01,yend = df.annot$ymin[5]) + annotate("segment", x=df.annot$xmax[5], xend=df.annot$xmax[5],y=df.annot$ymin[5]-0.01,yend = df.annot$ymin[5]) +
  annotate("segment",x=df.annot$xmin[6], xend=df.annot$xmin[6], y=df.annot$ymin[6]-0.01,yend = df.annot$ymin[6]) + annotate("segment", x=df.annot$xmax[6], xend=df.annot$xmax[6],y=df.annot$ymin[6]-0.01,yend = df.annot$ymin[6]) +
  annotate("segment",x=df.annot$xmin[7], xend=df.annot$xmin[7], y=df.annot$ymin[7]-0.01,yend = df.annot$ymin[7]) + annotate("segment", x=df.annot$xmax[7], xend=df.annot$xmax[7],y=df.annot$ymin[7]-0.01,yend = df.annot$ymin[7]) +
  annotate("segment",x=df.annot$xmin[8], xend=df.annot$xmin[8], y=df.annot$ymin[8]-0.01,yend = df.annot$ymin[8]) + annotate("segment", x=df.annot$xmax[8], xend=df.annot$xmax[8],y=df.annot$ymin[8]-0.01,yend = df.annot$ymin[8]) +
  annotate("segment",x=df.annot$xmin[9], xend=df.annot$xmin[9], y=df.annot$ymin[9]-0.01,yend = df.annot$ymin[9]) + annotate("segment", x=df.annot$xmax[9], xend=df.annot$xmax[9],y=df.annot$ymin[9]-0.01,yend = df.annot$ymin[9]) +
  annotate("segment",x=df.annot$xmin[10], xend=df.annot$xmin[10], y=df.annot$ymin[10]-0.01,yend = df.annot$ymin[10]) + annotate("segment", x=df.annot$xmax[10], xend=df.annot$xmax[10],y=df.annot$ymin[10]-0.01,yend = df.annot$ymin[10]) +
  annotate("segment",x=df.annot$xmin[11], xend=df.annot$xmin[11], y=df.annot$ymin[11]-0.01,yend = df.annot$ymin[11]) + annotate("segment", x=df.annot$xmax[11], xend=df.annot$xmax[11],y=df.annot$ymin[11]-0.01,yend = df.annot$ymin[11]) +
  annotate("segment",x=df.annot$xmin[12], xend=df.annot$xmin[12], y=df.annot$ymin[12]-0.01,yend = df.annot$ymin[12]) + annotate("segment", x=df.annot$xmax[12], xend=df.annot$xmax[12],y=df.annot$ymin[12]-0.01,yend = df.annot$ymin[12]) +
  annotate("segment",x=df.annot$xmin[13], xend=df.annot$xmin[13], y=df.annot$ymin[13]-0.01,yend = df.annot$ymin[13]) + annotate("segment", x=df.annot$xmax[13], xend=df.annot$xmax[13],y=df.annot$ymin[13]-0.01,yend = df.annot$ymin[13]) +
  annotate("segment",x=df.annot$xmin[14], xend=df.annot$xmin[14], y=df.annot$ymin[14]-0.01,yend = df.annot$ymin[14]) + annotate("segment", x=df.annot$xmax[14], xend=df.annot$xmax[14],y=df.annot$ymin[14]-0.01,yend = df.annot$ymin[14]) +
  annotate("segment",x=df.annot$xmin[15], xend=df.annot$xmin[15], y=df.annot$ymin[15]-0.01,yend = df.annot$ymin[15]) + annotate("segment", x=df.annot$xmax[15], xend=df.annot$xmax[15],y=df.annot$ymin[15]-0.01,yend = df.annot$ymin[15]) +
  annotate("text", x = df.annot$x, y = df.annot$y, label = df.annot$signif)

# Cross-group analysis
crossgroup <- tmp[tmp$Population1!=tmp$Population2,]

p <- ggplot(crossgroup, aes(x=Population1,y=w.UniFrac)) + theme_bw() +
  geom_violin(fill="grey80", col=F) + geom_boxplot(width=0.3,outlier.shape = 21)+
  facet_wrap(~Population2,ncol = 2,scales = "free_x")+
  xlab("Population") + ylab("Weighted UniFrac Distance")+ggtitle("Inter-Population Distance")

print(p)

crossgroup$pair <- paste(crossgroup$Population1, crossgroup$Population2,sep = "---")

pairwise.wilcox.test(x=crossgroup$w.UniFrac, g = crossgroup$pair,paired = F, p.adjust.method = "fdr") -> pairwg
pairwg <- data.frame(pairwg$p.value)
pairwg$Pair1 <- rownames(pairwg)
pairwg <- melt(pairwg, id.vars = "Pair1",value.name = "pvalue",variable.name = "Pair2")
pairwg <- na.omit(pairwg)
pairwg$Pair2 <- gsub(pattern = "\\.",replacement = "-",ignore.case = T,x = pairwg$Pair2)

pairwg$signif <- ifelse(pairwg$pvalue <0.001, "***",
                       ifelse(pairwg$pvalue <0.01, "**",
                              ifelse(pairwg$pvalue <0.05, "*",
                                     "ns")))

pairwg <- pairwg[order(pairwg$Pair1),]
print(pairwg)

balipairwg <- pairwg[grep(pattern = "Balinese---",x = pairwg$Pair1),]
balipairwg <- balipairwg[grep(pattern = "Balinese---",x = balipairwg$Pair2),]
balipairwg$Pair1 <- gsub(pattern = "Balinese---",replacement = "",ignore.case = T,x = balipairwg$Pair1)
balipairwg$Pair2 <- gsub(pattern = "Balinese---",replacement = "",ignore.case = T,x = balipairwg$Pair2)
balipairwg <- balipairwg[order(balipairwg$Pair1),]
rownames(balipairwg) <- NULL
balipairwg

df.annot <- balipairwg
rownames(df.annot) <- NULL
df.annot$Pair1 <- factor(df.annot$Pair1, levels = c("Hadza","Malawian","Guahibo","Italian","American"))
df.annot$Pair2 <- factor(df.annot$Pair2, levels = c("Hadza","Malawian","Guahibo","Italian","American"))

df.annot <- df.annot[order(df.annot$Pair2),]
df.annot <- df.annot[order(df.annot$Pair1),]


increase <- 0.025
df.annot$x <- (as.numeric(df.annot$Pair1) + as.numeric(df.annot$Pair2))/2
df.annot$y <- 0.3+seq(from=0.05, to=(increase*(1+nrow(df.annot))), by=increase)

df.annot$ymin <- df.annot$y-0.005
df.annot$ymax <- df.annot$y-0.005

df.annot$xmin <- as.numeric(df.annot$Pair1)
df.annot$xmax <- as.numeric(df.annot$Pair2)

print(df.annot)

p <- ggplot(subset(crossgroup, crossgroup$Population2 =="Balinese"), aes(x=Population1,y=w.UniFrac)) + theme_bw() +
  geom_violin(fill="grey80", col=F) + geom_boxplot(width=0.3,outlier.shape = 21)+
  xlab("Population") + ylab("Weighted UniFrac Distance")+ggtitle("Distance to Balinese")

p +
  annotate("segment",x = df.annot$xmin[1], xend = df.annot$xmax[1], y = df.annot$ymin[1], yend=df.annot$ymax[1])+
  annotate("segment",x = df.annot$xmin[2], xend = df.annot$xmax[2], y = df.annot$ymin[2], yend=df.annot$ymax[2])+
  annotate("segment",x = df.annot$xmin[3], xend = df.annot$xmax[3], y = df.annot$ymin[3], yend=df.annot$ymax[3])+
  annotate("segment",x = df.annot$xmin[4], xend = df.annot$xmax[4], y = df.annot$ymin[4], yend=df.annot$ymax[4])+
  annotate("segment",x = df.annot$xmin[5], xend = df.annot$xmax[5], y = df.annot$ymin[5], yend=df.annot$ymax[5])+
  annotate("segment",x = df.annot$xmin[6], xend = df.annot$xmax[6], y = df.annot$ymin[6], yend=df.annot$ymax[6])+
  annotate("segment",x = df.annot$xmin[7], xend = df.annot$xmax[7], y = df.annot$ymin[7], yend=df.annot$ymax[7])+
  annotate("segment",x = df.annot$xmin[8], xend = df.annot$xmax[8], y = df.annot$ymin[8], yend=df.annot$ymax[8])+
  annotate("segment",x = df.annot$xmin[9], xend = df.annot$xmax[9], y = df.annot$ymin[9], yend=df.annot$ymax[9])+
  annotate("segment",x = df.annot$xmin[10], xend = df.annot$xmax[10], y = df.annot$ymin[10], yend=df.annot$ymax[10])+
  annotate("segment",x=df.annot$xmin[1], xend=df.annot$xmin[1], y=df.annot$ymin[1]-0.01,yend = df.annot$ymin[1]) + annotate("segment", x=df.annot$xmax[1], xend=df.annot$xmax[1], y=df.annot$ymin[1]-0.01,yend = df.annot$ymin[1]) +
  annotate("segment",x=df.annot$xmin[2], xend=df.annot$xmin[2], y=df.annot$ymin[2]-0.01,yend = df.annot$ymin[2]) + annotate("segment", x=df.annot$xmax[2], xend=df.annot$xmax[2], y=df.annot$ymin[2]-0.01,yend = df.annot$ymin[2]) +
  annotate("segment",x=df.annot$xmin[3], xend=df.annot$xmin[3], y=df.annot$ymin[3]-0.01,yend = df.annot$ymin[3]) + annotate("segment", x=df.annot$xmax[3], xend=df.annot$xmax[3], y=df.annot$ymin[3]-0.01,yend = df.annot$ymin[3]) +
  annotate("segment",x=df.annot$xmin[4], xend=df.annot$xmin[4], y=df.annot$ymin[4]-0.01,yend = df.annot$ymin[4]) + annotate("segment", x=df.annot$xmax[4], xend=df.annot$xmax[4],y=df.annot$ymin[4]-0.01,yend = df.annot$ymin[4]) +
  annotate("segment",x=df.annot$xmin[5], xend=df.annot$xmin[5], y=df.annot$ymin[5]-0.01,yend = df.annot$ymin[5]) + annotate("segment", x=df.annot$xmax[5], xend=df.annot$xmax[5],y=df.annot$ymin[5]-0.01,yend = df.annot$ymin[5]) +
  annotate("segment",x=df.annot$xmin[6], xend=df.annot$xmin[6], y=df.annot$ymin[6]-0.01,yend = df.annot$ymin[6]) + annotate("segment", x=df.annot$xmax[6], xend=df.annot$xmax[6],y=df.annot$ymin[6]-0.01,yend = df.annot$ymin[6]) +
  annotate("segment",x=df.annot$xmin[7], xend=df.annot$xmin[7], y=df.annot$ymin[7]-0.01,yend = df.annot$ymin[7]) + annotate("segment", x=df.annot$xmax[7], xend=df.annot$xmax[7],y=df.annot$ymin[7]-0.01,yend = df.annot$ymin[7]) +
  annotate("segment",x=df.annot$xmin[8], xend=df.annot$xmin[8], y=df.annot$ymin[8]-0.01,yend = df.annot$ymin[8]) + annotate("segment", x=df.annot$xmax[8], xend=df.annot$xmax[8],y=df.annot$ymin[8]-0.01,yend = df.annot$ymin[8]) +
  annotate("segment",x=df.annot$xmin[9], xend=df.annot$xmin[9], y=df.annot$ymin[9]-0.01,yend = df.annot$ymin[9]) + annotate("segment", x=df.annot$xmax[9], xend=df.annot$xmax[9],y=df.annot$ymin[9]-0.01,yend = df.annot$ymin[9]) +
  annotate("segment",x=df.annot$xmin[10], xend=df.annot$xmin[10], y=df.annot$ymin[10]-0.01,yend = df.annot$ymin[10]) + annotate("segment", x=df.annot$xmax[10], xend=df.annot$xmax[10],y=df.annot$ymin[10]-0.01,yend = df.annot$ymin[10]) +
  annotate("text", x = df.annot$x, y = df.annot$y, label = df.annot$signif)

