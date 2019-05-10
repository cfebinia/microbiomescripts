rm(list=ls())

library(pvclust)
library(psych)

# load data
L6_CAG <- read.csv("bali_L6_common_taxa_transposed.csv",header = T,row.names = 1)

# Order by Abundance
n <- colSums(L6_CAG)
n <- sort(n, decreasing = T)
n <- names(n)

L6_CAG <- L6_CAG[,n]

# Load BMI data
df <- read.csv("data/bali_antropometry.csv",header = T,row.names = 1)
df <- df[rownames(L6_CAG),]
df$obesity <- ifelse(df$BMI >= 30, "obese","lean")

# Hierarchial Clustering on Pairwise Correlation Data on Lean Dataset

# calculate pairwise genus correlation
calc.rcorr <- corr.test(L6_CAG,method = "kendall",adjust = "fdr")

L6.kendall.cor <- round(calc.rcorr$r,5)
#write.table(L6.kendall.cor, file = "bali_genus_correlation_matrix_kendall.txt",quote = F,row.names = T)
#write.table(calc.rcorr$p, file = "bali_genus_correlation_p.values_matrix_kendall.txt",quote = F,row.names = T)

kendall.dist <- 1 - L6.kendall.cor

# Hierarchial Clustering Analysis (Ward.D)
hierclust <- hclust(as.dist(kendall.dist),method = "ward.D2")
plot(hierclust, hang = -1,
     main="Ward D2 clustering based on correlation distance (1 - Kendall) on Genus Abundance",cex=0.8)
rect.hclust(hierclust, k=10) -> clustk7

hclust.genus.order <- hierclust$labels[hierclust$order]
CAG.cat <- cutree(hierclust,k = 10)
CAG.cat <- CAG.cat[hclust.genus.order]

library(vegan)
cag.temp <- kendall.dist[hclust.genus.order,hclust.genus.order]

# Test Significance with ANOSIM
anoCAG <- anosim(cag.temp,grouping = CAG.cat,permutations = 10000)
anoCAG

sink("Figures/ANOSIM_CAG_K=10.txt")
print(anoCAG)
print(summary(anoCAG))
sink()

# Obtain distance
temp <- cmdscale(as.dist(cag.temp))
plot(temp)
coord.plot <- temp

# Obtain Dominant Taxon in each CAG
CAG.df.list <- list(CAG1=L6_CAG[,names(CAG.cat)[CAG.cat==1]],
                    CAG2=L6_CAG[,names(CAG.cat)[CAG.cat==2]],
                    CAG3=L6_CAG[,names(CAG.cat)[CAG.cat==3]],
                    CAG4=L6_CAG[,names(CAG.cat)[CAG.cat==4]],
                    CAG5=L6_CAG[,names(CAG.cat)[CAG.cat==5]],
                    CAG6=L6_CAG[,names(CAG.cat)[CAG.cat==6]],
                    CAG7=L6_CAG[,names(CAG.cat)[CAG.cat==7]],
                    CAG8=L6_CAG[,names(CAG.cat)[CAG.cat==8]],
                    CAG9=L6_CAG[,names(CAG.cat)[CAG.cat==9]],
                    CAG10=L6_CAG[,names(CAG.cat)[CAG.cat==10]])

cag.order <- names(CAG.df.list)

cag.dom <- unlist(lapply(CAG.df.list, function(x) names(x)[order(colSums(x),decreasing = T)][1]))
cag.dom <- as.character(cag.dom)
cag.dom

# Write CAG RA data by Genus
for(i in names(CAG.df.list)){
  write.csv(data.frame(CAG.df.list[[i]]),file = paste("data/bali_coabundance_RA_",i,".csv",sep=""), quote = F,col.names = F,row.names = T)
}

# Write CAG RA Data
CAG.RA <- lapply(CAG.df.list, function(x) rowSums(x))
CAG.RA <- data.frame(CAG.RA)
write.table(CAG.RA, "data/bali_coabundance_RA_byCAG.csv",quote = F,row.names = T,sep = ",")


# Abundance Heatmap
RAmap <- CAG.RA[rownames(df),]
RAmap <- data.frame(SampleID=rownames(RAmap),obesity=df$obesity, RAmap)

unifrac.clust <- read.csv("bali_cluster_ward's_rare90k.csv")
unifrac.clust$Order <- NULL
unifrac.clust$Cluster <- factor(unifrac.clust$Cluster, c("Type-P","Type-B"))

RAmap <- merge.data.frame(RAmap, unifrac.clust, by="SampleID")

apply(RAmap[,c(3:12)], 2, function(x) wilcox.test(x ~ RAmap$Cluster, paired = F,exact = F)) -> cagwilcox
do.call(rbind,lapply(cagwilcox, function(x) unlist(x[c("statistic","p.value")]))) -> cagwilcox

cagwilcox <- data.frame(cagwilcox)
cagwilcox$p.value <- round(cagwilcox$p.value, 3)
cagwilcox$signif <- ifelse(cagwilcox$p.value < 0.001, "***",
                           ifelse(cagwilcox$p.value < 0.01, "**",
                                  ifelse(cagwilcox$p.value < 0.05, "*",
                                         "ns")))

library(reshape2)
RAmap <- melt(RAmap, id.vars = c("SampleID","obesity", "Cluster"),variable.name = "CAG",value.name = "RA")
head(RAmap)

library(ggplot2)
ggplot(RAmap, aes(x=SampleID, y=CAG)) + theme_classic() +
  geom_tile(aes(fill = RA),lwd=0) +
  #scale_fill_gradient2(low = "dodgerblue3",mid = "aliceblue",high = "red3",midpoint = 10, limits=c(,70)) +
  scale_fill_gradient2(low = "dodgerblue",mid = "white",high = "red3",midpoint = 20,limits=c(0,60)) +
  facet_grid(~Cluster+obesity,space = "free",scales = "free_x") +
  theme(axis.text.x = element_text(angle=90,vjust = 0.3,hjust = 1,size=10,colour = "black")) +
  theme(axis.text.y = element_text(size=10,colour = "black"))+
  theme(axis.line = element_line(size=1.2,colour = "black"))+
  theme(axis.ticks = element_line(size=1.2,colour = "black"))+
  ylab("") + xlab("")+
  theme(legend.position = "bottom") +
  ggtitle("Fig3b -  CAG Abundance by Samples")

library(dplyr)

# std.err function
sem <- function(x){
  sd(x)/sqrt(length(x))
}

RAmap %>% group_by(CAG,Cluster) %>% summarise(RAmean=mean(RA), stderr=sem(RA)) -> ggbardata

p <- ggplot(ggbardata, aes(x=CAG, y=RAmean, fill=Cluster)) + theme_bw() +
  geom_bar(stat="identity", position=position_dodge(width = 0.5), width=0.5, col="black") +
  geom_errorbar(aes(group=Cluster,ymin=RAmean-2*stderr, ymax=RAmean+2*stderr), width=0.2, position=position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("grey40","grey80")) +
  scale_y_continuous(expand = c(0,0),limits = c(0,50))+
  theme(axis.text.x = element_text(size=10,colour = "black")) +
  theme(axis.text.y = element_text(size=10,colour = "black"))+
  theme(panel.grid = element_blank())+
  ylab("") + xlab("")+
  theme(legend.position = "right") +
  ggtitle("Fig3c - Differences in CAG Abundance between Community Types")

print(cagwilcox)
cagwilcox$x <- 1:10
cagwilcox$xmin <- cagwilcox$x -0.15
cagwilcox$xmax <- cagwilcox$x +0.15

cagwilcox$y <- c(35,40,20,20,15,15,30,15,20,15)

cagwilcox$ymin <- cagwilcox$y -1
cagwilcox$ymax <- cagwilcox$y

p + 
  annotate("segment",x=cagwilcox$xmin, xend=cagwilcox$xmax, y=cagwilcox$y, yend=cagwilcox$y) +
  annotate("segment",x=cagwilcox$xmin, xend=cagwilcox$xmin, y=cagwilcox$ymin, yend=cagwilcox$ymax) +
  annotate("segment",x=cagwilcox$xmax, xend=cagwilcox$xmax, y=cagwilcox$ymin, yend=cagwilcox$ymax) +
  annotate("text",label=cagwilcox$signif, x=cagwilcox$x, y=cagwilcox$y+1.5)



# memberships
library(dplyr)
library(reshape2)

CAGmemRA <- data.frame(CAG=c(), Taxa=c(), 
                       All_n=c(), All_MeanCI=c(),
                       TypeP_n=c(), TypeP_MeanCI=c(),
                       TypeB_n=c(), TypeB_MeanCI=c())

for(i in 1:10){
  x <- CAG.df.list[[i]]
  x[x==0] <- NA
  x$SampleID <- rownames(x)
  x <- merge.data.frame(x, unifrac.clust, by="SampleID")
  x <- melt(x, id.vars = c("SampleID","Cluster"),variable.name = "Taxa",value.name = "RA")
  x <- na.omit(x)
  x %>% group_by(Taxa) %>% summarise(All_n=length(RA),All_MeanCI=paste(round(mean(RA),2)," (",
                                                                       round((mean(RA)-1.96*(sd(RA)/sqrt(length(RA)))),2), " - ",
                                                                       round((mean(RA)+1.96*(sd(RA)/sqrt(length(RA)))),2), ")",
                                                                       sep="")) -> all
  x %>% group_by(Taxa, Cluster) %>% summarise(n=length(RA),MeanCI=paste(round(mean(RA),2)," (",
                                                                        round((mean(RA)-1.96*(sd(RA)/sqrt(length(RA)))),2), " - ",
                                                                        round((mean(RA)+1.96*(sd(RA)/sqrt(length(RA)))),2), ")",
                                                                        sep="")) -> bytypes
  bytypes$Cluster <- as.character(bytypes$Cluster)
  RA_typeP <- subset(bytypes, Cluster=="Type-P")
  colnames(RA_typeP)[c(3,4)] <- c("TypeP_n","TypeP_MeanCI")
  RA_typeB <- subset(bytypes, Cluster=="Type-B")
  colnames(RA_typeB)[c(3,4)] <- c("TypeB_n","TypeB_MeanCI")
  
  temp <- data.frame(CAG=paste("CAG",i, sep=""), all, RA_typeP[,c(3,4)], RA_typeB[,c(3,4)])
  
  CAGmemRA <- rbind(CAGmemRA, temp)
}

CAGmemRA

notCAG <- 100-rowSums(L6_CAG)
summary(notCAG)
paste(round(mean(notCAG),2)," (",
      round((mean(notCAG)-1.96*(sd(notCAG)/sqrt(length(notCAG)))),2), " - ",
      round((mean(notCAG)+1.96*(sd(notCAG)/sqrt(length(notCAG)))),2), ")",
      sep="") -> Other.all # Mean(95%CI intervals) 

n <- which(names(notCAG) %in% unifrac.clust$SampleID[unifrac.clust$Cluster=="Type-P"])
notCAG_P <- notCAG[n]
paste(round(mean(notCAG[n]),2)," (",
      round((mean(notCAG[n])-1.96*(sd(notCAG[n])/sqrt(length(notCAG[n])))),2), " - ",
      round((mean(notCAG[n])+1.96*(sd(notCAG[n])/sqrt(length(notCAG[n])))),2), ")",
      sep="") -> Other.P # Mean(95%CI intervals) for Type-P

n <- which(names(notCAG) %in% unifrac.clust$SampleID[unifrac.clust$Cluster=="Type-B"])
notCAG_B <- notCAG[n]
paste(round(mean(notCAG[n]),2)," (",
      round((mean(notCAG[n])-1.96*(sd(notCAG[n])/sqrt(length(notCAG[n])))),2), " - ",
      round((mean(notCAG[n])+1.96*(sd(notCAG[n])/sqrt(length(notCAG[n])))),2), ")",
      sep="") -> Other.B # Mean(95%CI intervals) for Type-B

wilcox.test(notCAG_P, notCAG_B,paired = F,exact = F) -> other.wilcox

Other <- data.frame(CAG="NA", Taxa="Other", All_n=41, All_MeanCI=Other.all, TypeP_n=41, TypeP_MeanCI=Other.P, 
                    TypeB_n=41, TypeB_MeanCI=Other.B, pval=as.character(round(other.wilcox$p.value,3)), fdr="ns")

# significant differences
common_compare <- L6_CAG
common_compare$SampleID <- rownames(common_compare)
common_compare <- merge.data.frame(unifrac.clust, common_compare, by="SampleID")
colnames(common_compare)

common_compare_wilcox <- apply(common_compare[,-c(1:2)], 2, function(x) wilcox.test(x ~ common_compare$Cluster, exact=F, paired=F))
common_pval <- lapply(common_compare_wilcox, function(x) x$p.value)
common_fdr <- p.adjust(unlist(common_pval), "fdr")
round(unlist(common_pval[(which(common_pval < 0.05))]),3)
names(which(common_fdr < 0.05))

CAGmemRA$pval <- 0
CAGmemRA$fdr <- 0

for(i in names(common_pval)){
  CAGmemRA[CAGmemRA$Taxa == i,"pval"] <- round(common_pval[[i]],3)
  CAGmemRA[CAGmemRA$Taxa == i,"fdr"] <- round(common_fdr[[i]],3)
}

CAGmemRA$pval <- ifelse(CAGmemRA$pval < 0.001, "0", CAGmemRA$pval)
CAGmemRA$fdr <- ifelse(CAGmemRA$fdr < 0.001, "0", CAGmemRA$fdr)

CAGmemRA$pval <- ifelse(CAGmemRA$pval >= 0.05, "ns", CAGmemRA$pval)
CAGmemRA$fdr <- ifelse(CAGmemRA$fdr >= 0.05, "ns", CAGmemRA$fdr)

CAGmemRA$pval <- ifelse(CAGmemRA$pval == "0", "<0.001", CAGmemRA$pval)
CAGmemRA$fdr <- ifelse(CAGmemRA$fdr == "0", "<0.001", CAGmemRA$fdr)

CAGmemRA <- rbind(CAGmemRA, Other)
CAGmemRA


write.table(CAGmemRA, "Figures/TableS6_bali_CAG_memberships&RA.csv", quote = F,row.names = T,sep = ",")


# heatmap
library(reshape2)
library(ggplot2)

df.cor <- L6.kendall.cor
df.cor <- df.cor[hclust.genus.order,rev(hclust.genus.order)]
df.cor <- as.data.frame(df.cor)
df.cor$Genus1 <- rownames(df.cor)
df.cor <- melt(df.cor, id.vars = "Genus1",variable.name = "Genus2",value.name = "kendall")
df.cor$Genus1 <- factor(df.cor$Genus1, levels = hclust.genus.order)
df.cor$Genus2 <- factor(df.cor$Genus2, levels = rev(hclust.genus.order))

p<- ggplot(df.cor, aes(x=Genus1, y=Genus2)) + theme_classic() +
  geom_tile(aes(fill = kendall),lwd=0) +
  scale_fill_gradient2(low = "#15399A",mid = "white",high = "red3",midpoint = 0, limit=c(-1,1)) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.3,hjust = 1,size=10,colour = "black")) +
  theme(axis.text.y = element_text(size=10,colour = "black"))+
  theme(axis.line = element_line(size=1.2,colour = "black"))+
  theme(axis.ticks = element_line(size=1.2,colour = "black"))+
  ylab("") + xlab("")+
  theme(legend.position = "top") +
  ggtitle("FigS3 - Kendall Correlationship Matrix")

print(p)

# Create Wiggum Plot (network)
library(GGally)

# define network color
nodecol <- c("blue","cyan2","gold","magenta","salmon","forestgreen","red","dodgerblue","chartreuse3","grey30")
names(nodecol) <- cag.order
nodecol

#select only significant and positive  relationship
select.cor <- data.frame(L6.kendall.cor[hclust.genus.order,hclust.genus.order])
select.cor$tax1 <- rownames(select.cor)
select.cor <- melt(select.cor,id.vars = "tax1",variable.name = "tax2", value.name = "Kendall_t")
head(select.cor)
pval <- data.frame(calc.rcorr$p[hclust.genus.order,hclust.genus.order])
pval$tax1 <- rownames(pval)
pval <- melt(pval,id.vars = "tax1",variable.name = "tax2", value.name="pval")
head(pval)
select.cor$pval <- pval$pval

head(select.cor)

select.cor$print <- "YES"
select.cor$print[select.cor$pval > 0.05] <- "NO"
select.cor$print[select.cor$value <= 0] <- "NO"
select.cor$print[which(select.cor$tax1 == select.cor$tax2)] <- "NO"

network.df <- subset(select.cor, print=="YES")
head(network.df)

#create network object
net <- network(network.df[,c(1,2)],directed = F)

ggnet2(net)

# Set color by CAG
Taxacols <- CAGmemRA[,c("CAG","Taxa")]
cols <- data.frame(CAG=names(nodecol), color=nodecol,stringsAsFactors = F)
Taxacols <- merge.data.frame(Taxacols, cols, by="CAG")
Taxacols <- na.omit(Taxacols)
Taxacols$CAG <- as.character(Taxacols$CAG)
Taxacols$Taxa <- as.character(Taxacols$Taxa)

nodes.data <- ggnet2(net,mode = "fruchtermanreingold")$data
nodes.data$CAG <- ""

for(i in Taxacols$Taxa){
  nodes.data[which(nodes.data$label ==i), "CAG"] <- Taxacols$CAG[Taxacols$Taxa==i]
  nodes.data[which(nodes.data$label ==i), "color"] <- Taxacols$color[Taxacols$Taxa==i]
}

# set group positions
GroupByVertex02 = function(Groups) {
  numGroups = length(unique(Groups))
  GAngle    = (1:numGroups) * 2 * pi / numGroups
  Centers   = matrix(c(cos(GAngle), sin(GAngle)), ncol=2)/3
  x = y = c()
  for(i in 1:numGroups) {
    curGroup = which(Groups == unique(Groups)[i])
    VAngle = (1:length(curGroup)) * 2 * pi / length(curGroup)
    x = c(x, Centers[i,1] + cos(VAngle) / (numGroups))
    y = c(y, Centers[i,2] + sin(VAngle) / (numGroups))
  }
  matrix(c(x, y), ncol=2)
}

Groups <- nodes.data$color
GBV2 <- GroupByVertex02(Groups)

names <- c()
for(i in 1:length(unique(Groups))){
  names <- c(names,nodes.data$label[nodes.data$color==unique(Groups)[i]])
}
names

rownames(GBV2) <- names
colnames(GBV2) <- c("x","y")

net %v% "x" = GBV2[nodes.data$label,"x"]
net %v% "y" = GBV2[nodes.data$label,"y"]

ggnet2(net, mode=c("x","y"),node.color = "nodcols", edge.color = c("color", "grey"),node.size = 3)

#set node size
nodes.data$size <- colSums(L6_CAG[,nodes.data$label]*90000)
nodes.data$size <- nodes.data$size/sum(nodes.data$size)*1000
nodes.data$size[nodes.data$size<=20] <- 2

node.size <- sqrt(nodes.data$size)*2
node.size[node.size < 10] <- 4

ggnet2(net, mode=c("x","y"),node.color = "nodcols", edge.color = c("color", "grey"),
       node.size = nodes.data$size) +
  theme(legend.position = "none")

# set font face
nodes.data$fontface = "plain"
nodes.data[which(nodes.data$label %in% cag.dom),"fontface"] <- "bold"
nodes.data$fontsize = 2
nodes.data[which(nodes.data$label %in% cag.dom),"fontsize"] <- 5

# set node positions
library(vegan)
t <- metaMDS(kendall.dist)
t <- as.matrix(t$points)
t <- t[nodes.data$label,]
net %v% "x.nmds" = t[,1]
net %v% "y.nmds" = t[,2]


#set node colors
customcols <- as.character(nodes.data$color)
set.vertex.attribute(net,"nodcols",value = customcols)

# draw network
ggnet2(net, mode=c("x.nmds","y.nmds"),color = "nodcols", palette = unique(nodes.data$color), size=0,
       edge.color = "grey", edge.size=1, layout.exp = 0.5) +
  geom_point(aes(color = color), size = node.size, color = "grey40") +
  geom_point(aes(color = color), size = node.size, alpha = 0.5) +
  geom_point(aes(color = color), size = node.size-2) +
  geom_text(aes(label = nodes.data$label, fontface=nodes.data$fontface),size=5, color = "black") +
  guides(color = FALSE)


# Display by CAG only
nodecol

for(i in nodecol){
  set.vertex.attribute(net,"nodcols",value = ifelse(customcols == i, i,"grey80"))
  textcol <- ifelse(customcols == i, "black",NA)
  
  # draw network
  p <- ggnet2(net, mode=c("x.nmds","y.nmds"),color = "nodcols", palette = unique(nodes.data$color), size=0,
              edge.color = "grey80", edge.size=1, layout.exp = 0.5) +
    geom_point(aes(color = color), size = node.size, color = "grey80") +
    geom_point(aes(color = color), size = node.size, alpha = 0.5) +
    geom_point(aes(color = color), size = node.size-2) +
    geom_text(aes(label = nodes.data$label, fontface=nodes.data$fontface), color = textcol, size=6) +
    guides(color = FALSE)
  
  print(p)
}

# Write CAG Abundance
L6_CAG <- t(L6_CAG)
L6_CAG <- data.frame(Taxon=rownames(L6_CAG),L6_CAG)
rownames(L6_CAG) <- NULL

colnames(Taxacols)[2] <- "Taxon"

L6_CAG <- merge.data.frame(Taxacols[,-3],L6_CAG, by="Taxon")
head(L6_CAG)

write.csv(L6_CAG,file = "data/bali_common_Genera_L6_by_CAG.csv",quote = F,col.names = F,row.names = F)
