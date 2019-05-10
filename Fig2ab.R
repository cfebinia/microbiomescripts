#### *Prevotella* and *Bacteroides* distribution across 6 populations

rm(list=ls())

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

plot(cap,type="n")
points(cap, disp="sites", cex=1.5, pch=c(15,2,1)[gam.data$Lifestyle])
legend("topright",pch=c(15,2,1),levels(gam.data$Lifestyle),xjust = 1,yjust=1, cex=1.1)

# load Relative abundance matrix
L6_taxa <- read.csv("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/6_population_L6_Relative_abundance_tidy.csv",header=T)  

L6RA_PB <- L6_taxa[L6_taxa$Genus %in% c("Prevotella","Bacteroides"),-c(1:5)]
rownames(L6RA_PB) <- L6RA_PB$Genus
L6RA_PB <- t(L6RA_PB[,-1])
L6RA_PB <- L6RA_PB[sample_name,]


#### Generalized Linear Model analysis on beta diversity data (Supplementary Results)
library(mgcv)
library(sp)
library(plotrix)

#create palette
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"),space="Lab",interpolate="linear")

#extract MDS axes
w.uni.pcoa.axes <- summary(cap)$sites

# create GAM dataset
cols <- c("SampleID","Country","Location","Lifestyle",
          "Age","Sex","BMI","Obesity","MDS1","MDS2","MDS3")


gam.data <- cbind(environ.var[rownames(w.uni.pcoa.axes),colnames(environ.var) %in% cols],
                  w.uni.pcoa.axes[,colnames(w.uni.pcoa.axes) %in% cols],
                  Prevotella=L6RA_PB[,"Prevotella"], Bacteroides = L6RA_PB[,"Bacteroides"])

# find convex in the PCoA plot
no.cols<-256
gr<-101				# resolution along each axis. 
findConvex<-function(x,y,names){
  hull<-cbind(x,y)[chull(cbind(x,y)),]
  px<-pretty(x)			# returns same number of items as x, but rounded and evenly spaced.
  py<-pretty(y)
  x.new<-seq(min(px),max(px),len=gr)
  y.new<-seq(min(py),max(py),len=gr)
  ingrid<-as.data.frame(expand.grid(x.new,y.new))
  Fgrid<-ingrid
  Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
  names(Fgrid)<-names
  return(Fgrid)
}

# Pre-processing. These are the coordinates at which GAM predictions are extracted, for plotting surfaces. Can't plot surfaces
dfgenus <- findConvex(gam.data$MDS1, gam.data$MDS2, c("MDS1","MDS2"))
dfgenus$MDS3 <- median(gam.data$MDS3)

##### Fig. 2a - *Prevotella* abundance fitted to weighted UniFrac data with GAM
name <- "Prevotella"

training <- gam.data
newdata <- dfgenus

k1 <- 10
k2 <- 10

gam1 <- gam(eval(parse(text=name)) ~ s(MDS1,k=k1,bs="tp")+s(MDS2,k=k1,bs="tp")+s(MDS3,k=k1,bs="tp") 
            +s(MDS1,MDS2,k=k2,bs="tp")+s(MDS1,MDS3,k=k2,bs="tp") +s(MDS2,MDS3,k=k2,bs="tp")
            +s(MDS1,MDS2,MDS3,k=k2,bs="tp"), 
            gamma=1.0,method="REML",data=training,select=TRUE)

fit1 <- predict(gam1, newdata=newdata)

nlev <- 8
map <- rgb.palette(no.cols)

# get the minimum and maximum values enountered in the fit. min(c(unlist(gam.data[11:12])), na.rm=T) #
mn <- min(c(unlist(fit1)),na.rm=TRUE)		
mx <- max(c(unlist(fit1)),na.rm=TRUE)

locs<-(range(unlist(fit1),na.rm=TRUE)-mn)/(mx-mn)*no.cols
surf<-matrix(fit1,nrow=sqrt(dim(newdata)[1]))
px<-pretty(gam.data$MDS1)			# the tick locations, made pretty (0  5 10 15 20 25 30)
py<-pretty(gam.data$MDS2)
x.new<-seq(min(px),max(px),len=gr)	# x and y locations where coloured boxes are to be drawn. 
y.new<-seq(min(py),max(py),len=gr)
xlab <- axis_lab[1]
ylab <- axis_lab[2]

par(mar=c(5,5,5,5))
image(x.new,y.new,surf,col=map[locs[1]:locs[2]],xlab=xlab,ylab=ylab,
      main=name,sub="PCoA of Weighted Unifrac Distance",axes=FALSE)
axis(1)								# draw the axes
axis(2)
abline(h=0,v=0,lty="dashed",col="grey40")
abline(h=1,v=1,col="black")
points(cap, disp="sites", cex=1.3, pch=c(15,2,1)[gam.data$Lifestyle])
legend("topright",pch=c(15,2,1),levels(gam.data$Lifestyle),xjust = 1,yjust=1, cex=1.1)

sink("Figures/Prevotella.GAM.txt")
print(name)
print(summary(gam1))
sink()

##### Fig. 2b - *Bacteroides* abundance fitted on weighted UniFrac data  with GAM
name <- "Bacteroides"

training <- gam.data
newdata <- dfgenus

k1 <- 10
k2 <- 10

gam1 <- gam(eval(parse(text=name)) ~ s(MDS1,k=k1,bs="tp")+s(MDS2,k=k1,bs="tp")+s(MDS3,k=k1,bs="tp") 
            +s(MDS1,MDS2,k=k2,bs="tp")+s(MDS1,MDS3,k=k2,bs="tp") +s(MDS2,MDS3,k=k2,bs="tp")
            +s(MDS1,MDS2,MDS3,k=k2,bs="tp"), 
            gamma=1.0,method="REML",data=training,select=TRUE)

fit1 <- predict(gam1, newdata=newdata)

nlev <- 8
map <- rgb.palette(no.cols)

# get the minimum and maximum values enountered in the fit. min(c(unlist(gam.data[11:12])), na.rm=T) #
mn <- min(c(unlist(fit1)),na.rm=TRUE)		
mx <- max(c(unlist(fit1)),na.rm=TRUE)

locs<-(range(unlist(fit1),na.rm=TRUE)-mn)/(mx-mn)*no.cols
surf<-matrix(fit1,nrow=sqrt(dim(newdata)[1]))
px<-pretty(gam.data$MDS1)			# the tick locations, made pretty (0  5 10 15 20 25 30)
py<-pretty(gam.data$MDS2)
x.new<-seq(min(px),max(px),len=gr)	# x and y locations where coloured boxes are to be drawn. 
y.new<-seq(min(py),max(py),len=gr)
xlab <- axis_lab[1]
ylab <- axis_lab[2]

par(mar=c(5,5,5,5))
image(x.new,y.new,surf,col=map[locs[1]:locs[2]],xlab=xlab,ylab=ylab, main=name,sub="PCoA of Weighted Unifrac Distance",axes=FALSE)
axis(1)								# draw the axes
axis(2)
abline(h=0,v=0,lty="dashed",col="grey40")
abline(h=1,v=1,col="black")
points(cap, disp="sites", cex=1.3, pch=c(15,2,1)[gam.data$Lifestyle])
legend("topright",pch=c(15,2,1),levels(gam.data$Lifestyle),xjust = 1,yjust=1, cex=1.1)

sink("Figures/Bacteroides.GAM.txt")
print(name)
print(summary(gam1))
sink()

##### Fig. 2c - Distribution of *Prevotella* and *Bacteroides* abundance across 6 populations

boxdata <- data.frame(Population=as.character(environ.var[sample_name,"Population"]), L6RA_PB[sample_name,])

boxdata <- melt(boxdata, id.vars = c("Population"),variable.name = "Genus",value.name = "RelativeAbundance")

boxdata$Population <- factor(boxdata$Population, c("Hadza","Malawian","Guahibo","Balinese","American","Italian"))
boxdata$Genus <- factor(boxdata$Genus, c("Prevotella","Bacteroides"))

boxlabel <- summary(boxdata$Population)/2

boxdata$xlab <- ""

for(i in names(boxlabel)){
  boxdata[boxdata$Population==i,"xlab"] <- paste(i,"\n(n=",boxlabel[i],")", sep="")
}

boxdata <- boxdata[order(boxdata$Population),]
boxdata$xlab <- factor(boxdata$xlab, levels = unique(boxdata$xlab))

ggplot(boxdata, aes(x=xlab, y=RelativeAbundance, fill=Genus))+ theme_bw() +
  geom_boxplot(col = "black",outlier.shape = 21, outlier.size = 2, width=0.5)+
  scale_fill_manual(values = c("grey80","grey40")) +
  ylab("Relative Abundance (%)") + xlab("") +
  ggtitle("Fig 2c - Distribution of Prevotella and Bacteroides Abundance Across Populations")
