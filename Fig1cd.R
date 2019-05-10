### Fig. 1c- Fig. 1d

rm(list=ls())

library(vegan)
library(stats)
library(dplyr)
library(kableExtra)

#load distance matrix
bali.w.unifrac.dm <- read.table("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/bali_weighted_unifrac_rarefaction_90000_0.txt",header = T)

#reorder samples
bali.sample.ord <- sort(rownames(bali.w.unifrac.dm))
bali.w.unifrac.dm <- bali.w.unifrac.dm[bali.sample.ord,bali.sample.ord]

# draw PCoA
bali.w.unifrac.dist <- as.dist(bali.w.unifrac.dm)
cap <- capscale(bali.w.unifrac.dist ~ 1)
eig <- cap$CA$eig
var <- eig/sum(eig)
axis_labs <- paste("PCo",1:length(names(var)), " (", round(var*100,1),"%)", sep="")
axis_labs

# Hierarchial clustering of weighted UniFrac data (Ward's method)
bali.cluster <- hclust(bali.w.unifrac.dist, method = "ward.D2")
plot(bali.cluster, hang=-1)
rect.hclust(bali.cluster, k=2)

# Load BMI data
df <- read.csv("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/bali_antropometry.csv",header = T,row.names = 1)
df <- df[bali.sample.ord,]
df$obesity <- ifelse(df$BMI >= 30, "obese","lean")

# Cluster picks, using k=2 to divide the tree in two groups 
cluster <- cutree(bali.cluster, k=2)
cluster1 <- names(which(cluster==1))
cluster2 <- names(which(cluster==2))

df$cluster <- "NA"
df[cluster1,]$cluster <- "Cluster 1"
df[cluster2,]$cluster <- "Cluster 2"

df$obesity <- factor(df$obesity)
df$cluster <- factor(df$cluster, c("Cluster 1","Cluster 2"))


# Apply the info into the previous ordination plot
plot(cap, type="n",xlab=axis_labs[1], ylab=axis_labs[2],cex.axis=1.7, cex.lab=1.7)
points(cap, display="sites", pch=c(21,24)[factor(df$obesity)], 
       #bg=c("navy","turquoise3")[df$cluster], 
       bg=NA,
       col=c("black","turquoise4")[df$cluster], 
       cex=c(1.8,2.2)[factor(df$obesity)])
title("Fig 1c - Weighted UniFrac with Ward Clusters")

# ANOSIM
types.anosim<- anosim(dat = bali.w.unifrac.dist,grouping =  df$cluster,permutations = 100000)
types.anosim


### Fig. 1d

#load taxonomy data
#L6_tax <- read.csv("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/bali_L6_commonGenera.csv",header = T,row.names = 1)

L6_tax <- read.csv("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/bali_L6_common_taxa_transposed.csv",header = T,row.names = 1)

# sort by abundance
L6_abund <- round(L6_tax*1e6)
reads.bySample <- rowSums(L6_tax)
reads.bygenus <- colSums(L6_tax) 
#all(((reads.bygenus*100/sum(L6_tax)) >= 0.01) == TRUE)  # check that all taxa > 1%

# sort SampleID
L6_tax <- L6_tax[bali.sample.ord,]

# Vector regression analysis
L6.envfit <- envfit(cap,L6_tax, permutations = 100000) 
L6.envfit.df <- data.frame(r2=L6.envfit$vectors$r,pval=L6.envfit$vectors$pvals)

L6.envfit.df$significance <- "ns"
L6.envfit.df[L6.envfit.df$pval < 0.0001,"significance"] <- "***"

L6.envfit.signif <- L6.envfit.df[L6.envfit.df$pval< 0.0001,]

# display result
#kable(L6.envfit.signif, caption="Table S2 Vector regression analysis of weighted UniFrac distance against bacterial genus abundance (only taxa with p <0.0001 are shown).", align = c("l",rep('c',6))) %>%
#  kable_styling(full_width = F)

L6.envfit.df

n <- rownames(L6.envfit.signif)
L6.permanova <- adonis(bali.w.unifrac.dist ~ ., L6_tax[,rownames(L6.envfit.signif)],permutations = 100000)
print(L6.permanova)

sink("Figures/Cluster_ANOSIM_VR_PERMANOVA.txt")
print(summary(types.anosim))
print(L6.envfit.df)
print(L6.permanova)
sink()

# Select the taxon of interest
constrain.taxa <- L6.permanova$aov.tab
constrain.taxa <- rownames(constrain.taxa)[constrain.taxa$`Pr(>F)` < 0.0001]
constrain.taxa <- constrain.taxa[which(constrain.taxa!="NA")]
constrain.taxa <- L6_tax[,constrain.taxa]

# Constrained ordination
ccap <- capscale(bali.w.unifrac.dist ~ ., data=constrain.taxa)
eig <- c(ccap$CCA$eig, ccap$CA$eig)
var <- eig/sum(eig)
axis_labs <- paste(names(var), " (", round(var*100,1),"%)", sep="")
axis_labs

#draw the plot
plot(ccap, type="n",xlab=axis_labs[1], ylab=axis_labs[2],cex.axis=1.7, cex.lab=1.7)
points(ccap, display="sites", 
       pch=c(21,24)[factor(df$obesity)], 
       #bg=c("navy","turquoise3")[df$cluster], 
       bg=NA,
       col=c("black","turquoise4")[df$cluster], 
       cex=c(1.8,2.2)[factor(df$obesity)])
text(ccap, display="bp",cex=1,col="blue")
title("Fig 1d - Constrained Weighted UniFrac")


# Write cluster Info
bali.cluster <- hclust(bali.w.unifrac.dist, method = "ward.D2")
sample.order.clust <- bali.cluster$labels[bali.cluster$order]

bali.cluster <- data.frame(SampleID=rownames(df), Cluster=as.character(df$cluster),stringsAsFactors = F)
bali.cluster$Cluster[bali.cluster$Cluster=="Cluster 1"] <- "Type-P"
bali.cluster$Cluster[bali.cluster$Cluster=="Cluster 2"] <- "Type-B"


bali.cluster$SampleID <- factor(bali.cluster$SampleID, levels = sample.order.clust)
bali.cluster$Order <- as.numeric(bali.cluster$SampleID)

write.csv(bali.cluster,file="data/bali_cluster_ward's_rare90k.csv",quote = F,col.names = F,row.names = F)

