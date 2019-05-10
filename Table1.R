### Description of studied populations (Table 1)
library(dplyr)
library(reshape2)

#Load environmental data
environ.var <- read.csv("C:/Users/Clarissa/Google Drive/FKMUI-EIJKMAN Nutriepid/data paper/Manuscript/v5/Additional Files/published dataset/6_populations_metadata_new.csv", header=TRUE, row.names=1,stringsAsFactors = T)

# setting factorial variables
popord <- c("Balinese","Hadza","Malawian","Guahibo","Italian","American")
environ.var$Population <- factor(environ.var$Population, popord)
environ.var$Obesity <- ifelse(is.na(environ.var$BMI), "NA", ifelse(environ.var$BMI > 30, "obese","lean"))
environ.var$Obesity <- factor(environ.var$Obesity, c("lean","obese","NA"))
environ.var$Lifestyle <- factor(environ.var$Lifestyle, c("Transitional","Industrial","PreIndustrial"))
environ.var$Sex <- factor(environ.var$Sex, c("Male","Female"))
environ.var$femcount <- as.numeric(environ.var$Sex) -1

# summarise data
environ.var %>% group_by(Population) %>% 
  summarise(Sample.Size=length(Age), 
            Age.mean=as.integer(mean(Age,na.rm = T)),
            Age.sd=as.integer(sd(Age,na.rm = T)),
            Female.Percentage=paste(as.integer(sum(femcount,na.rm = T)*100/length(femcount)),"%",sep=""),
            Sequence.Count.Mean=as.integer(mean(SequenceCount)),
            Sequencing.Region=paste(unique(SequencingRegion)),
            Platform=paste(unique(Platform))) -> table1

print(table1)
write.csv(table1, file="Figures/Table1.csv",row.names = F, col.names = F,quote = F)
