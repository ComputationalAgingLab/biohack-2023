library(readr)
metareg_37 <- read_csv("Expression/metareg_37.csv")


library(GEOquery)

#gsm <- getGEO(filename="GSE40279/GSE40279_family.soft.gz")
gsm <- getGEO(filename="GSE557663/GSE55763_family.soft.gz")
Illumina450k<-gsm@gpls$GPL13534@dataTable@table
#Illumina450k$ID[grep(metareg_37$...1[1],Illumina450k$UCSC_RefGene_Name)]

NNl<-Illumina450k[grep(metareg_37$...1[1],Illumina450k$UCSC_RefGene_Name),][1,]
NNl[1,]<-NA
Anatation_Gene_CpG<-list()
for(i in 1:nrow(metareg_37)){
  
  A<-Illumina450k[grep(metareg_37$...1[i],Illumina450k$UCSC_RefGene_Name),]
  if(nrow(A) == 0){A<-NNl}
  
  Anatation_Gene_CpG[[i]]<-cbind(metareg_37$...1[i],A)
  
  
}
Anatation_Gene_CpG2<-do.call(rbind,Anatation_Gene_CpG)

row.names(top.table1)[top.table1$adj.P.Val<0.05]

Anatation_Gene_CpG2_1<-Anatation_Gene_CpG2[Anatation_Gene_CpG2$ID %in% row.names(top.table1)[top.table1$adj.P.Val<0.05], ]

Anatation_Gene_CpG2_2<-Anatation_Gene_CpG2[Anatation_Gene_CpG2$ID %in% row.names(top.table2)[top.table2$adj.P.Val<0.05], ]

Anatation_Gene_CpG2_3<-Anatation_Gene_CpG2[Anatation_Gene_CpG2$ID %in% row.names(top.table3)[top.table3$adj.P.Val<0.05], ]

load("C:/Hakaton/Aging/GSE557663_GSE40279_GSE152027_beta.RData")

sampleTable$gender<-gsub(" ","",sampleTable$gender)


sampleTable$age<-round(sampleTable$age)

sampleTable$group<-2
sampleTable$group[sampleTable$age<37]<-1
sampleTable$group[sampleTable$age>60]<-3
sampleTable$group<-as.factor(sampleTable$group)

library(ggplot2)

CpG2_1<-GSE557663_GSE40279_GSE152027_beta[Anatation_Gene_CpG2_1$ID,]
CpG2_2<-GSE557663_GSE40279_GSE152027_beta[Anatation_Gene_CpG2_2$ID,]

CpG2_3<-GSE557663_GSE40279_GSE152027_beta[unique(c(Anatation_Gene_CpG2_1$ID,Anatation_Gene_CpG2_2$ID)),]


save(CpG2_1,CpG2_2,Anatation_Gene_CpG2_1,Anatation_Gene_CpG2_2,sampleTable,file = "Reduct_CpG.RData")


for( i in 1:nrow(CpG2_1)){
  tpm<-data.frame(beta = t(CpG2_1[i,]),sampleTable)
  
  
  p <- ggplot(tpm, aes(age,tpm[,1], color=group, shape=group))
  png(file=paste0("pic_G_1/",colnames(tpm)[1],"_",Anatation_Gene_CpG2_1[Anatation_Gene_CpG2_1$ID == colnames(tpm)[1],1],".png"), width=4096, height=2048, res=300)
  
  plot(p + geom_point() + geom_smooth(method=lm, aes(fill=group))+ theme_classic())
  dev.off();
}

tpm<-data.frame(M = log2(t(CpG2_1[1,])/(1-t(CpG2_1[1,]))),sampleTable)


p <- ggplot(tpm, aes(tpm[,1], age, color=group, shape=group))
p + geom_point() + geom_smooth(method=lm, aes(fill=group))+ theme_classic()




for( i in 1:nrow(CpG2_2)){
  tpm<-data.frame(beta = t(CpG2_2[i,]),sampleTable)
  
  
  p <- ggplot(tpm, aes(age,tpm[,1], color=group, shape=group))
  png(file=paste0("pic_G_2/",colnames(tpm)[1],"_",Anatation_Gene_CpG2_2[Anatation_Gene_CpG2_2$ID == colnames(tpm)[1],1],".png"), width=4096, height=2048, res=300)
  
  plot(p + geom_point() + geom_smooth(method=lm, aes(fill=group))+ theme_classic())
  dev.off();
}



write.table(data.frame(Anatation_Gene_CpG2_1$ID),row.names = F,col.names = F, quote = F,file = "Anatation_Gene_CpG2_1.txt")

write.table(data.frame(Anatation_Gene_CpG2_1$`metareg_37$...1[i]`),row.names = F,col.names = F, quote = F,file = "Anatation_Gene_CpG2_1_Gene.txt")


write.table(data.frame(Anatation_Gene_CpG2_2$ID),row.names = F,col.names = F, quote = F,file = "Anatation_Gene_CpG2_2.txt")

write.table(data.frame(Anatation_Gene_CpG2_2$`metareg_37$...1[i]`),row.names = F,col.names = F, quote = F,file = "Anatation_Gene_CpG2_2_Gene.txt")

