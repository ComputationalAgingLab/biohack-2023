library(EpiDISH)
library(GEOquery)


gsm <- getGEO(filename="GSE40279/GSE40279_family.soft.gz")
gsm@header
names(gsm@gsms)
GSE40279_Sample<-do.call(rbind,lapply(gsm@gsms, function(x){unlist(lapply(strsplit(x@header$characteristics_ch1,split = ":"),function(x){x[2]}))}))
colnames(GSE40279_Sample)<-c("age","source","plate","gender","ethnicity","tissue")
GSE40279_Sample<-as.data.frame(GSE40279_Sample)
GSE40279_Sample$age<-as.numeric(GSE40279_Sample$age)


library(readr)
GSE40279_average_beta_txt <- read_delim("C:/Hakaton/GSE40279_average_beta.txt.gz", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

GSE40279_average_beta_txt<-as.data.frame(GSE40279_average_beta_txt)
row.names(GSE40279_average_beta_txt)<-GSE40279_average_beta_txt$ID_REF
GSE40279_average_beta_txt<-GSE40279_average_beta_txt[,-1]
colnames(GSE40279_average_beta_txt)<-row.names(GSE40279_Sample)

save(GSE40279_Sample,GSE40279_average_beta_txt,file = "GSE40279.RData")

################################################



gsm <- getGEO(filename="GSE557663/GSE55763_family.soft.gz")

GSE557663_Sample<-do.call(rbind,lapply(gsm@gsms, function(x){unlist(lapply(strsplit(x@header$characteristics_ch1,split = ":"),function(x){x[2]}))}))
gsm@gsms
GSE557663_Sample$name2<-unlist(lapply(gsm@gsms, function(x){unlist(lapply(strsplit(x@header$title,split = ", "),function(x){x[2]}))}))

GSE557663_Sample2<-GSE557663_Sample[match(colnames(GSE557663_average_beta_txt),GSE557663_Sample$name2),]
GSE557663_Sample<-GSE557663_Sample2
sum(GSE557663_Sample2$name2==colnames(GSE557663_average_beta_txt))
colnames(GSE557663_average_beta_txt)[1]
colnames(GSE557663_Sample)<-c("tissue","dataset","gender","age")

GSE557663_average_beta_txt <- read_delim("C:/Hakaton/GSE55763_normalized_betas.txt.gz", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)
COLNa<-GSE557663_average_beta_txt[,1]
Beta_filtr<-GSE557663_average_beta_txt[,seq(3,ncol(GSE557663_average_beta_txt),2)] > 0.05
GSE557663_average_beta_txt<-GSE557663_average_beta_txt[,seq(2,ncol(GSE557663_average_beta_txt),2)]
GSE557663_average_beta_txt<-as.data.frame(GSE557663_average_beta_txt)
summary(colSums(Beta_filtr))
row.names(GSE557663_average_beta_txt)<-COLNa$ID_REF
quantile(colSums(Beta_filtr),probs = 0.998)
quantile(rowSums(Beta_filtr),probs = 0.9999)
sum(Beta_filtr)/(ncol(Beta_filtr)*nrow(Beta_filtr))
sum(rowSums(Beta_filtr) < (ncol(Beta_filtr)*0.05))
GSE557663_average_beta_txt<-GSE557663_average_beta_txt[(rowSums(Beta_filtr) < (ncol(Beta_filtr)*0.05)),]



save(GSE557663_Sample,GSE557663_average_beta_txt,file = "GSE557663.RData")

save(GSE557663_Sample,file = "GSE557663_Sample.RData")


#########################################################################

gsm <- getGEO(filename="GSE152027_family.soft.gz")
gsm@gsms

GSE152027_Sample<-do.call(rbind,lapply(gsm@gsms, function(x){unlist(lapply(strsplit(x@header$characteristics_ch1,split = ": "),function(x){x[2]}))}))
colnames(GSE152027_Sample)<-c("status","gender","age")
GSE152027_Sample<-as.data.frame(GSE152027_Sample)
GSE152027_Sample$age<-as.numeric(GSE152027_Sample$age)

GSE152027_Sample$group<-2
GSE152027_Sample$group[GSE152027_Sample$age<37]<-1
GSE152027_Sample$group[GSE152027_Sample$age>60]<-3
GSE152027_Sample$name2<-unlist(lapply(gsm@gsms, function(x){unlist(lapply(strsplit(x@header$title,split = " "),function(x){x[1]}))}))

GSE152027_Sample$name2 == colnames(GSE152027_signals)




library(readr)
GSE152027_signals <- read_csv("GSE152027_IOP_processed_signals.csv.gz")
GSE152027_signals<-as.data.frame(GSE152027_signals)
row.names(GSE152027_signals)<-GSE152027_signals$...1
GSE152027_signals<-GSE152027_signals[,-1]
GSE152027_signals<-GSE152027_signals[,seq(1,ncol(GSE152027_signals),2)]


GSE152027_Sample2<-GSE152027_Sample[match(colnames(GSE152027_signals),GSE152027_Sample$name2),]
sum(GSE152027_Sample$name2 == colnames(GSE152027_signals))

GSE557663_Sample<-GSE152027_Sample2

frac.m_GSE152027<- epidish(beta.m = GSE152027_signals, ref.m = centDHSbloodDMC.m,method = 'RPC')


save(GSE152027_Sample,frac.m_GSE152027,GSE152027_signals,file = "GSE152027.RData")


#####################################################################################


load("/mnt/SHD_Users/Zarubin/Aging/GSE557663.RData")
load("/mnt/SHD_Users/Zarubin/Aging/GSE557663_Sample.RData")
load("/mnt/SHD_Users/Zarubin/Aging/GSE557663_age_cells.RData")


load("/mnt/SHD_Users/Zarubin/Aging/GSE40279.RData")
load("/mnt/SHD_Users/Zarubin/Aging/GSE40279_age_cells.RData")



GSE40279_Sample$group<-2
GSE40279_Sample$group[GSE40279_Sample$age<37]<-1
GSE40279_Sample$group[GSE40279_Sample$age>60]<-3

sampleTable1<-cbind(GSE557663_Sample[,c(4,3,6)],as.data.frame(frac.m_GSE557663$estF))
sampleTable2<-cbind(GSE40279_Sample[,c(1,4,7)],as.data.frame(frac.m_GSE40279$estF))
sampleTable1$Exp<-"GSE557663"
sampleTable2$Exp<-"GSE40279"

sampleTable<-rbind(sampleTable1,sampleTable2)

sampleTable$group<-as.factor(sampleTable$group)

GSE557663_average_beta_txt<-GSE557663_average_beta_txt[row.names(GSE557663_average_beta_txt) %in% row.names(GSE40279_average_beta_txt),]
GSE40279_average_beta_txt<-GSE40279_average_beta_txt[row.names(GSE40279_average_beta_txt) %in% row.names(GSE557663_average_beta_txt),]
GSE40279_average_beta_txt<-GSE40279_average_beta_txt[match(row.names(GSE557663_average_beta_txt),row.names(GSE40279_average_beta_txt)),]
sum(row.names(GSE557663_average_beta_txt) == row.names(GSE40279_average_beta_txt))


GSE557663_GSE40279_beta<-cbind(GSE557663_average_beta_txt,GSE40279_average_beta_txt)
save(sampleTable,GSE557663_GSE40279_beta,file = "GSE557663_GSE40279_beta.RData")

GSE152027_Sample<-GSE152027_Sample[GSE152027_Sample$status %in% "CON",]
GSE152027_Sample<-GSE152027_Sample[!is.na(GSE152027_Sample$age),]
GSE152027_Sample$group<-2
GSE152027_Sample$group[GSE152027_Sample$age<37]<-1
GSE152027_Sample$group[GSE152027_Sample$age>60]<-3

#as.data.frame(frac.m_GSE152027$estF)[GSE152027_Sample$name2,]
sampleTable3<-cbind(GSE152027_Sample[,c(3,2,4)],as.data.frame(frac.m_GSE152027$estF)[GSE152027_Sample$name2,])

GSE152027_signals<-GSE152027_signals[,GSE152027_Sample$name2]
colnames(GSE152027_signals)<-row.names(GSE152027_Sample)
GSE152027_signals<-GSE152027_signals[row.names(GSE152027_signals) %in% row.names(GSE557663_GSE40279_beta),]
GSE557663_GSE40279_beta<-GSE557663_GSE40279_beta[row.names(GSE557663_GSE40279_beta) %in% row.names(GSE152027_signals),]
sum(row.names(GSE152027_signals) == row.names(GSE557663_GSE40279_beta))
GSE557663_GSE40279_beta<-GSE557663_GSE40279_beta[match(row.names(GSE152027_signals),row.names(GSE557663_GSE40279_beta)),]
sum(row.names(GSE152027_signals) == row.names(GSE557663_GSE40279_beta))

sampleTable3$Exp<-"GSE152027"
sampleTable<-rbind(sampleTable,sampleTable3)
GSE557663_GSE40279_GSE152027_beta<-cbind(GSE557663_GSE40279_beta,GSE152027_signals)
sampleTable$gender<-gsub(" ","",sampleTable$gender)
sampleTable$group<-as.factor(sampleTable$group)
save(sampleTable,GSE557663_GSE40279_GSE152027_beta,file = "GSE557663_GSE40279_GSE152027_beta.RData")



