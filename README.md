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
write.csv(resdata, file="differ_1.csv")
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

## WGCNA
tf_expr <- read.csv("tf_expression.csv", row.names = 1, header = TRUE, check.names = FALSE)
pathway_expr <- read.csv("pathway_expression.csv", row.names = 1, header = TRUE, check.names = FALSE)


common_samples <- intersect(colnames(tf_expr), colnames(pathway_expr))
tf_expr <- tf_expr[, common_samples]
pathway_expr <- pathway_expr[, common_samples]

datExpr <- t(rbind(tf_expr, pathway_expr))


if(!is.numeric(datExpr[1, 1])) {
  datExpr <- apply(datExpr, 2, as.numeric)
}
library(WGCNA)
options(stringsAsFactors = FALSE)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes] # 移除不好的基因和样本
}

sampleTree <- hclust(dist(datExpr), method = "average")
par(mar = c(3, 3, 2, 1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.90, col="red") 
net <- blockwiseModules(datExpr,
                        power = 6, 
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendrogram = FALSE,
                        saveTOMs = FALSE,
                        verbose = 3)

table(net$colors)
moduleColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pathway_genes_in_datExpr <- colnames(datExpr)[colnames(datExpr) %in% rownames(pathway_expr)]
pathway_expression <- rowMeans(datExpr[, pathway_genes_in_datExpr, drop = FALSE])
traitData <- data.frame(PathwayExpression = pathway_expression)
rownames(traitData) <- rownames(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0) 
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(3, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Pathway Membership",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-Pathway Relationships")
module <- "turquoise"
moduleGenes <- (moduleColors == module) 
cat("Number of genes in turquoise module:", sum(moduleGenes), "\n")
cat("Length of moduleMembership:", length(moduleMembership), "\n")
cat("Length of GeneSignificance:", length(GeneSignificance), "\n")
cat("Length of moduleMembership[moduleGenes]:", length(moduleMembership[moduleGenes]), "\n")
cat("Length of GeneSignificance[moduleGenes]:", length(GeneSignificance[moduleGenes]), "\n")

adjacency <- adjacency(datExpr, power = 6)
intramodularConnectivity <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)
modConnectivity <- intramodularConnectivity[moduleGenes, "kWithin"]
modConnectivity <- sort(modConnectivity, decreasing = TRUE)
head(modConnectivity, 10)
hub_genes <- names(head(modConnectivity, 20))
tf_hubs <- hub_genes[hub_genes %in% rownames(tf_expr)]
print(tf_hubs)
GS1 <- as.numeric(cor(datExpr, traitData$PathwayExpression, use = "p"))
GeneSignificance <- abs(GS1)
datKME <- signedKME(datExpr, MEs)
MM <- datKME

module <- "turquoise"
column <- match(paste0("kME", module), names(MM))
moduleMembership <- MM[, column]

par(mfrow = c(1,1))
verboseScatterplot(abs(moduleMembership[moduleGenes]),
                   GeneSignificance[moduleGenes],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for pathway",
                   main = paste("MM vs. GS\n"),
                   col = module)

turquoise_genes <- colnames(datExpr)[moduleColors == "turquoise"]
cat("Turquoise 模块中的基因数量:", length(turquoise_genes), "\n")
print(turquoise_genes)

## Z-score
library(readxl) 
library(pheatmap) 
library(RColorBrewer) 

data <- read_excel("MYB_lingin_count.xlsx")
data_df <- as.data.frame(data)
rownames(data_df) <- data_df[,1]
data_df <- data_df[,-1]
data_matrix <- as.matrix(data_df)
data_log <- log2(data_matrix + 1)
data_zscore <- t(scale(t(data_log)))

pheatmap(data_zscore,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         border_color = NA,
         main = "Gene Expression Heatmap (Z-score)")
zscore_df <- as.data.frame(data_zscore)
zscore_df$Gene <- rownames(zscore_df) 
zscore_df <- zscore_df[, c(ncol(zscore_df), 1:(ncol(zscore_df)-1))]
write.csv(zscore_df, "MYB_lingin_zscore_data.csv", row.names = FALSE)

