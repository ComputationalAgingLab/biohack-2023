library(lumi)
library(limma)



design<- model.matrix(~  0 + group +   .,data = sampleTable)
#colnames(design)[4]<-"gender_M"

fit <- lmFit(GSE557663_GSE40279_GSE152027_beta, design)

tmp <- contrasts.fit(fit, coefficients = "age")
tmp <- eBayes(tmp)

top.table_age <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_age$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "genderM")
tmp <- eBayes(tmp)

top.table_gender <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_age$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "B")
tmp <- eBayes(tmp)

top.table_B <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_B$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "NK")
tmp <- eBayes(tmp)

top.table_NK <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_NK$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "CD4T")
tmp <- eBayes(tmp)

top.table_CD4T <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_CD4T$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "CD8T")
tmp <- eBayes(tmp)

top.table_CD8T <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_CD8T$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "Mono")
tmp <- eBayes(tmp)

top.table_Mono <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_Mono$adj.P.Val<0.05)


tmp <- contrasts.fit(fit, coefficients = "Neutro")
tmp <- eBayes(tmp)

top.table_Neutro <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_Neutro$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "ExpGSE40279")
tmp <- eBayes(tmp)

top.table_ExpGSE40279 <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_ExpGSE40279$adj.P.Val<0.05)

tmp <- contrasts.fit(fit, coefficients = "ExpGSE557663")
tmp <- eBayes(tmp)

top.table_ExpGSE557663 <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table_ExpGSE557663$adj.P.Val<0.05)


contr <- makeContrasts(group1 - group2, levels = colnames(design))

tmp <- contrasts.fit(fit, contr)

#tmp <- contrasts.fit(fit, coefficients = "group1")

#Функция eBayes применяется для статистического тестирования на основе усиленного метода Байеса.

tmp <- eBayes(tmp)

top.table1 <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table1$adj.P.Val<0.05)

contr <- makeContrasts(group1 - group3, levels = colnames(design))

tmp <- contrasts.fit(fit, contr)
#Функция eBayes применяется для статистического тестирования на основе усиленного метода Байеса.

tmp <- eBayes(tmp)

top.table2 <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table2$adj.P.Val<0.05)


contr <- makeContrasts(group2 - group3, levels = colnames(design))

tmp <- contrasts.fit(fit, contr)

#Функция eBayes применяется для статистического тестирования на основе усиленного метода Байеса.

tmp <- eBayes(tmp)

top.table3 <- topTable(tmp, sort.by = "P", n = Inf)
sum(top.table3$adj.P.Val<0.05)


dif_group_cpg<-unique(c(row.names(top.table1)[top.table1$adj.P.Val<0.05],
                        row.names(top.table2)[top.table2$adj.P.Val<0.05],
                        row.names(top.table3)[top.table3$adj.P.Val<0.05]))
length(dif_group_cpg)

dif_group_cpg_filt<-dif_group_cpg[!(dif_group_cpg %in%  row.names(top.table_gender)[top.table_gender$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_B)[top.table_B$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_CD4T)[top.table_CD4T$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_CD8T)[top.table_CD8T$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_Mono)[top.table_Mono$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_Neutro)[top.table_Neutro$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_NK)[top.table_NK$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_ExpGSE40279)[top.table_ExpGSE40279$adj.P.Val<0.05])]
length(dif_group_cpg_filt)
dif_group_cpg_filt<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_ExpGSE557663)[top.table_ExpGSE557663$adj.P.Val<0.05])]
length(dif_group_cpg_filt)

dif_group_cpg_filt2<-dif_group_cpg_filt[!(dif_group_cpg_filt %in%  row.names(top.table_age)[top.table_age$adj.P.Val<0.05])]
length(dif_group_cpg_filt2)

save(fit,top.table_age,top.table1,top.table2,top.table3,dif_group_cpg_filt,file = "top.table_B.RData")


library(pheatmap)

df <- sampleTable[,1:3]
row.names(df)<-colnames(GSE557663_GSE40279_GSE152027_beta)



Tmp<-GSE557663_GSE40279_GSE152027_beta[dif_group_cpg_filt,]
plot(t(Tmp)[,1],sampleTable$age)

#dB<-rowMeans(Tmp[,sampleTable$group == 1],na.rm = T) - rowMeans(Tmp[,sampleTable$group == 2],na.rm = T)
#sum(abs(dB)>0.05)
#Tmp2<-Tmp[abs(dB)>0.05,]

pheatmap(Tmp, cluster_rows=T, show_rownames=F,show_colnames=F,
         cluster_cols=T,annotation_col  = df)

df2 <- sampleTable[,1:3]
ORe<-order(df2$age)
Tmp2<-GSE557663_GSE40279_GSE152027_beta[dif_group_cpg_filt,ORe]
df2 <- sampleTable[ORe,1:3]
row.names(df2)<-colnames(Tmp2)

pheatmap(Tmp2, cluster_rows=T, show_rownames=F,show_colnames=F,
         cluster_cols=F,annotation_col  = df2)