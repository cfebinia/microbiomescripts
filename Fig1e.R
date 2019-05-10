# Taxa Family Bar Charts (Fig 1e)
#setwd("C:/Users/Clarissa/Google Drive/Bali_diet_paper/")
library(reshape2)
library(dplyr)
library(ggplot2)

#load L5 otu data
df_L5 <- read.table(file="bali_L5_taxa_summary_relativeabund_parsed_rare90k.txt",header =T, stringsAsFactors = F)
head(df_L5)
dim(df_L5)

# Tidy Unclassified
L5_ID <- df_L5$Family
n <- which(L5_ID == "Unclassified")
L5_ID[n] <- paste(df_L5$Order,"OTU", sep=".")[n]
df_L5$Family <- L5_ID

# Combine Other with Low Abundance Taxa
df_L5 <- df_L5[,-(which(colnames(df_L5)%in% c("Kingdom","Class","Order")))]

L5_Other <- df_L5$Family
n <- which(L5_Other == "Other")
L5_Other <- df_L5[n,]
dim(L5_Other)

df_L5 <- df_L5[-n,]
L5_LowAbund <- rowSums(df_L5[-c(1,2)])
names(L5_LowAbund) <- df_L5$Family
n <- names(which(L5_LowAbund < 1))
n <- which(df_L5$Family %in% n)

L5_LowAbund <- df_L5[n ,]
dim(L5_LowAbund)
df_L5 <- df_L5[-n,]

L5_Other <- rbind(L5_Other, L5_LowAbund)
dim(L5_Other)

L5_Other <- cbind(data.frame(Phylum="Bacteria.Other",Family="Other"),
                  t(data.frame(colSums(L5_Other[,-c(1:2)]))))

rownames(L5_Other) <- NULL

df_L5 <- rbind(L5_Other,df_L5)
dim(df_L5)

df_L5$Family <- as.character(df_L5$Family)
df_L5$Phylum <- as.character(df_L5$Phylum)

# Melt and Sort by Abundance
temp <-melt(df_L5, id.vars = c("Phylum","Family"),variable.name = "SampleID",value.name = "RelativeAbundance",factorsAsStrings = F)
temp %>% group_by(Family) %>% summarise(RA = sum(RelativeAbundance)) -> family.ord
family.ord <- family.ord[order(family.ord$RA,decreasing = T),]$Family
n <- which(family.ord %in% c("Prevotellaceae","Bacteroidaceae"))
family.ord[n] <- rev(family.ord[n])

temp %>% group_by(Phylum) %>% summarise(RA = sum(RelativeAbundance)) -> phylum.ord
phylum.ord <- phylum.ord[order(phylum.ord$RA,decreasing = T),]$Phylum

df_L5$Family <- factor(df_L5$Family, levels = family.ord)
df_L5$Phylum <- factor(df_L5$Phylum, levels = phylum.ord)

df_L5 <- df_L5[order(df_L5$Family),]
df_L5 <- df_L5[order(df_L5$Phylum),]
head(df_L5)

df_L5$Family <- as.character(df_L5$Family)
df_L5$Phylum <- as.character(df_L5$Phylum)

n <- grep(pattern = "Other",x =df_L5$Phylum)
df_L5[n,"Phylum"] <- "Bacteria.Other"

n <- grep(pattern = "Other",x =df_L5$Family)
df_L5 <- rbind(df_L5[-n,],df_L5[n,])

#set colors
mycols_by_phylum <- read.csv("color_setup_by_phylum.csv",stringsAsFactors = F)
mycols_by_phylum <- rbind(mycols_by_phylum, 
                          data.frame(Phylum=c("Firmicutes2","Bacteroidetes2","Proteobacteria2"), 
                                     cols=c("hotpink","royalblue4","darkolivegreen4")))
mycols_by_phylum[mycols_by_phylum$Phylum=="Bacteroidetes","cols"] <- "skyblue1"
mycols_by_phylum[mycols_by_phylum$Phylum=="Firmicutes","cols"] <- "orangered1"
mycols_by_phylum[mycols_by_phylum$Phylum=="Proteobacteria","cols"] <- "green3"
mycols_by_phylum[mycols_by_phylum$Phylum=="Verrucomicrobia","cols"] <- "white"
mycols_by_phylum[mycols_by_phylum$Phylum=="Bacteria.Other","cols"] <- "grey30"

famlist_Firmicutes <- subset(df_L5, df_L5$Phylum=="Firmicutes")$Family
famlist_Bacteroidetes<- subset(df_L5, df_L5$Phylum=="Bacteroidetes")$Family[-c(1:4)]
famlist_Proteobacteria<- subset(df_L5, df_L5$Phylum=="Proteobacteria")$Family

is.even <- function(x) x %% 2 == 0

n <- length(famlist_Firmicutes)
even_F <- famlist_Firmicutes[is.even(1:n)]

n <- length(famlist_Bacteroidetes)
even_B <- famlist_Bacteroidetes[is.even(1:n)]

n <- length(famlist_Proteobacteria)
even_P <- famlist_Proteobacteria[is.even(1:n)]

temp <- df_L5
temp$Phylum2 <- as.character(temp$Phylum)
temp[which(temp$Family %in% even_F),"Phylum2"] <- "Firmicutes2"
temp[which(temp$Family %in% even_B),"Phylum2"] <- "Bacteroidetes2"
temp[which(temp$Family %in% even_P),"Phylum2"] <- "Proteobacteria2"


mycolsF <- c()
for(i in unique(temp$Phylum2)){
  print(i)
  famlist <- temp[temp$Phylum2 == i,"Family"]
  n <- length(famlist)
  base <- mycols_by_phylum[mycols_by_phylum$Phylum ==i,"cols"]
  colfunc <- colorRampPalette(c(base, "#FFFFFF"))
  colist <- c(colfunc(n+1)[1:n])
  names(colist) <- famlist
  mycolsF <- c(mycolsF,colist)
}
mycolsF

mycolsF["Bacteroidaceae"] <- "blue"
mycolsF["Prevotellaceae"] <- "turquoise3"
mycolsF["Paraprevotellaceae*"] <- "turquoise1"
mycolsF["Rikenellaceae"] <- "steelblue1"
mycolsF["Veillonellaceae"] <-"firebrick1"
mycolsF["Ruminococcaceae"] <-"red4"
mycolsF["Lachnospiraceae"] <- "red2"
mycolsF["Enterobacteriaceae"] <- "forestgreen"

#load antrop data
#df_antrop <- read.table("antropometry.txt",header = T,stringsAsFactors = F)
#SampleID <- colnames(df_L5)[-c(1,2)]
#df_antrop <- df_antrop[which(df_antrop$SampleID %in% SampleID),]
#df_antrop$obesity <- ifelse(df_antrop$BMI >= 30, "obese","lean")

# make plot dataset
ggdata <- melt(df_L5, id.vars = c("Phylum","Family"),variable.name = "SampleID",value.name = "RelativeAbundance",factorsAsStrings = F)

# reorder by Phylum
phylum.ord <- unique(ggdata$Phylum)
phylum.ord <- phylum.ord[c(1,3:14,2)]
ggdata$Phylum <- factor(ggdata$Phylum, levels = phylum.ord)
ggdata <- ggdata[order(ggdata$Phylum),]

# reorder by Phylum
family.ord <- unique(ggdata$Family)
n <- which(family.ord %in% c("Bacteroidaceae","Prevotellaceae"))
family.ord <- c(family.ord[-n], c("Bacteroidaceae","Prevotellaceae"))

ggdata$Family <- factor(ggdata$Family, levels = family.ord)
ggdata <- ggdata[order(ggdata$Family),]

# reoder by sampleID cluster
bali.cluster <- read.csv("bali_cluster_ward's_rare90k.csv",stringsAsFactors = F)
sample.order.clust <- bali.cluster$Order
names(sample.order.clust) <- bali.cluster$SampleID
sample.order.clust <- sort(sample.order.clust)
sample.order.clust <- names(sample.order.clust)
ggdata$SampleID <- factor(ggdata$SampleID, levels = sample.order.clust)

ggplot(ggdata, aes(x=SampleID, y=RelativeAbundance, fill=Family)) +
         theme_minimal() + 
         geom_bar(stat="identity", width = 0.9) + 
         scale_fill_manual(values=mycolsF) +
         scale_y_continuous(expand = c(0,0)) +
         scale_x_discrete(expand = c(0,0)) +
         theme(axis.text.x = element_text(angle=90,vjust = 0.5)) +
         theme(axis.title.y = element_text(face = "bold")) +
         theme(legend.position = "bottom") +
         theme(legend.text=element_text(size=10)) +
         theme(panel.grid=element_blank()) +
         theme(axis.ticks.y = element_blank()) +
         theme(axis.ticks.x = element_blank()) +
         theme(axis.line.x = element_blank()) +
         theme(axis.line.y = element_blank()) +
         theme(panel.background = element_rect(fill="white",colour = "white")) +
         xlab("") + ylab("") +
  guides(fill= guide_legend(reverse = T))



