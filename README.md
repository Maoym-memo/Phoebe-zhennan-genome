# Phoebe-zhennan-genome

./OrthoFinder/orthofinder -f proteinS/longestpep/ -S diamond -M msa -A mafft -t 8 -a 8



## deseq2-pheatmap
setwd("D:/data/04 RNA-seq expression/new")
us_count<-read.csv("read.csv",head=T,row.names=1) 
us_count[is.na(us_count)] <- 0 
us_count<-round(us_count,digits=0) 

us_count<-as.matrix(us_count) 
condition<-factor(c("STEM","STEM","STEM","LEAF","LEAF","LEAF","CAMBIA","CAMBIA","CAMBIA","FLO","FLO","FLO","DRUPE","DRUPE","DRUPE"))
coldata<-data.frame(row.names=colnames(us_count),condition)  
coldata 
condition 
library(DESeq2) 
dds<-DESeqDataSetFromMatrix(us_count,coldata,design=~condition)
head(dds)
dds<-DESeq(dds) 
resultsNames(dds)  
res<-results(dds) 
summary(res)
plotMA(res,ylim=c(-2,2)) 
mcols(res,use.names=TRUE)
plot(res$log2FoldChange,res$pvalue) 
res <- res[order(res$padj),]
resdata <-merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
deseq_res<-data.frame(resdata)
up_diff_result<-subset(deseq_res,padj < 0.05 & (log2FoldChange > 1)) 
down_diff_result<-subset(deseq_res,padj < 0.05 & (log2FoldChange < -1)) 
write.csv(up_diff_result,"UPgene.csv") 
write.csv(down_diff_result,"DOWNgene.csv") 
write.csv(resdata, file="differ_1.csv")

#pheatmap绘制热图
library("pheatmap")
data <- read.csv("differ_1.csv", header = TRUE, row.names = 1)
data_matrix <- data[, c("Pzhe661", "Pzhe665", "Pzhe655R","Pzhe667", "Pzhe669", "Pzhe663", "Pzhe670", 
                        "Pzhe662", "Pzhe666", "Pzhe668", "Pzhe657", "Pzhe6611", "Pzhe660", 
                        "Pzhe664", "Pzhe659")]
head(data_matrix, n = 5)
data_filtered <- data_matrix[rowSums(data_matrix) > 0, ]
dim(data_filtered)
pheatmap(data_filtered, 
         scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         fontsize_row = 3, 
         fontsize_col = 10)
