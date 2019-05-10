#### Alpha Diversity by population (Fig. 1a)

library(dplyr)
library(reshape2)
library(ggplot2)


# std.err function
sem <- function(x){
  sd(x)/sqrt(length(x))
}

# Faith's PD
adiv.faith <- read.csv("6_population_alphadiv_PD_rarefaction_4000_new.csv",header = T,stringsAsFactors = F)

# summarise data
adiv.faith %>% group_by(sequences.per.sample,Population,Metric) %>% 
  summarise(Sample.Size=length(value), 
            Mean=round(mean(value),2), 
            StDev=round(sd(value),2),
            StErr=round(sem(value),2),
            LowerCI=round(mean(value)-2*sem(value),2),
            UpperCI=round(mean(value)+2*sem(value),2)) -> ggdata

# specific colors and order were assigned to each population.
ggdata$Legend <- paste(ggdata$Population, "\n(n=", ggdata$Sample.Size,")", sep="")
ggdata$Population <- factor(ggdata$Population, levels=c("Hadza","Malawian","Guahibo","Balinese","Italian","American"))
ggdata <- ggdata[order(ggdata$Population),]
ggdata$Legend <- factor(ggdata$Legend, levels=unique(ggdata$Legend))

popcols <- c("gold","mediumseagreen","forestgreen","red","skyblue","dodgerblue3")
names(popcols) <- unique(ggdata$Legend)

# draw plot
library(ggplot2, quietly = F,warn.conflicts = F)

ggplot(ggdata, aes(x=sequences.per.sample, y=Mean, col=Legend)) +
  theme_bw() + 
  geom_line(size=1) + theme(panel.grid.minor = element_blank())+ #theme(legend.position = c(0.85, 0.2)) +
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, col=Legend), size=0.5, width=100) +
  geom_point(shape=21, fill="white", size=3) +
  scale_color_manual(values = popcols) +
  scale_x_continuous(limits = c(0,4500), breaks = c(100,500,seq(1000,4500,by = 500)))+
  xlab("Sequences/Sample") + ylab("Faith's PD")+
  ggtitle("Phylogenetic Diversity")

# Simpson
adiv.simpson <- read.csv("6_population_alphadiv_simpson_rarefaction_4000.csv",header = T,stringsAsFactors = F)

# summarise data
adiv.simpson %>% group_by(sequences.per.sample,Population,Metric) %>% 
  summarise(Sample.Size=length(value), 
            Mean=round(mean(value),2), 
            StDev=round(sd(value),2),
            StErr=round(sem(value),2),
            LowerCI=round(mean(value)-2*sem(value),2),
            UpperCI=round(mean(value)+2*sem(value),2)) -> ggdata

# specific colors and order were assigned to each population.
ggdata$Legend <- paste(ggdata$Population, "\n(n=", ggdata$Sample.Size,")", sep="")
ggdata$Population <- factor(ggdata$Population, levels=c("Hadza","Malawian","Guahibo","Balinese","Italian","American"))
ggdata <- ggdata[order(ggdata$Population),]
ggdata$Legend <- factor(ggdata$Legend, levels=unique(ggdata$Legend))

popcols <- c("gold","mediumseagreen","forestgreen","red","skyblue","dodgerblue3")
names(popcols) <- unique(ggdata$Legend)

# draw plot
library(ggplot2, quietly = F,warn.conflicts = F)

ggplot(ggdata, aes(x=sequences.per.sample, y=Mean, col=Legend)) +
  theme_bw() + 
  geom_line(size=1) + theme(panel.grid.minor = element_blank())+ #theme(legend.position = c(0.85, 0.2)) +
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, col=Legend), size=0.5, width=100) +
  geom_point(shape=21, fill="white", size=3) +
  scale_color_manual(values = popcols) +
  scale_x_continuous(limits = c(0,4500), breaks = c(100,500,seq(1000,4500,by = 500)))+
  xlab("Sequences/Sample") + ylab("Simpson 1/D Index")+
  ggtitle("Alpha Diversity (Simpson)")

# Shannon
adiv.shannon <- read.csv("6_population_alphadiv_shannon_rarefaction_4000.csv",header = T,stringsAsFactors = F)

# summarise data
adiv.shannon %>% group_by(sequences.per.sample,Population,Metric) %>% 
  summarise(Sample.Size=length(value), 
            Mean=round(mean(value),2), 
            StDev=round(sd(value),2),
            StErr=round(sem(value),2),
            LowerCI=round(mean(value)-2*sem(value),2),
            UpperCI=round(mean(value)+2*sem(value),2)) -> ggdata

# specific colors and order were assigned to each population.
ggdata$Legend <- paste(ggdata$Population, "\n(n=", ggdata$Sample.Size,")", sep="")
ggdata$Population <- factor(ggdata$Population, levels=c("Hadza","Malawian","Guahibo","Balinese","Italian","American"))
ggdata <- ggdata[order(ggdata$Population),]
ggdata$Legend <- factor(ggdata$Legend, levels=unique(ggdata$Legend))

popcols <- c("gold","mediumseagreen","forestgreen","red","skyblue","dodgerblue3")
names(popcols) <- unique(ggdata$Legend)

# draw plot
library(ggplot2, quietly = F,warn.conflicts = F)

ggplot(ggdata, aes(x=sequences.per.sample, y=Mean, col=Legend)) +
  theme_bw() + 
  geom_line(size=1) + theme(panel.grid.minor = element_blank())+ #theme(legend.position = c(0.85, 0.2)) +
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI, col=Legend), size=0.5, width=100) +
  geom_point(shape=21, fill="white", size=3) +
  scale_color_manual(values = popcols) +
  scale_x_continuous(limits = c(0,4500), breaks = c(100,500,seq(1000,4500,by = 500)))+
  xlab("Sequences/Sample") + ylab("Shannon Index")+
  ggtitle("Alpha Diversity (Shannon)")
