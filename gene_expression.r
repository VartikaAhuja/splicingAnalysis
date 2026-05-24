## DESeq2 to shortlist upregulated splicing factor genes

## Load DESeq2 and input data
# Read gene count matrix (remove chromosome/start/end columns)
library(DESeq2)
library("RColorBrewer")
library("pheatmap")
getwd()
count_data=read.table("gene_counts.txt",sep="\t",header=TRUE,row.names="gene")
#the first 3 columns for gene count data have chromosome, start and end information for genes, when generated using bedtools
count_data = count_data[,-c(1,2,3)]
head(count_data)
# Read sample metadata (includes sample name, type, condition, etc.)
col_data=read.table(file = "samples_info.txt", header = T, sep = "\t")
rownames(col_data) = col_data$Sample

all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))

## Create DESeq2 dataset
dds = DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

# Filter out low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# Set 'normal' as reference condition
dds$condition <- relevel(dds$condition, ref = "normal")

## Transformation and visualization
rld = rlogTransformation(dds)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition))
colnames(mat) = NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# Generate heatmap and PCA plots for sample clustering
pdf(file="sample_pheatmap.pdf")
pheatmap(mat, clustering_distance_rows=distsRL, clustering_distance_cols=distsRL, col=colors)
dev.off()
pdf(file="sample_pca.pdf")
plotPCA(rld, intgroup=c("condition"))
dev.off()

## Differential expression analysis
dds = estimateSizeFactors(dds)
sizeFactors(dds)
dds = DESeq(dds)

res<-results(dds)
head(res)
resOrdered <- res[order(rownames(res)),]

# Save ordered and filtered results
write.table(resOrdered,file="deseq2_all_results.txt",sep="\t",quote=F,row.names=T, col.names=T)
resFiltered <- subset(resOrdered,resOrdered$pvalue<0.05)
write.table(resFiltered,file="deseq2_all_pval0.05.txt",sep="\t",quote=F,row.names=T, col.names=T)
#To check the upregulated genes among a pre-defined set of genes for study; where genelist has those genes list
dat<- read.table("genelist",sep="\t",header=F,col.names="gene")
resFiltered2 <- subset(resOrdered,rownames(resOrdered) %in% dat$gene)
dim(resFiltered)
dim(resFiltered2)
resFiltered3 <- subset(resFiltered2,resFiltered2$pvalue<0.05)
dim(resFiltered3)
pdf(file="plot_ma.pdf")
plotMA(res, ylim=c(-7,7))
dev.off()

# Identify upregulated/downregulated genes with logFC thresholds
upRes2=resFiltered3[resFiltered3$log2FoldChange>=2,]
downRes2=resFiltered3[resFiltered3$log2FoldChange<=-2,]
upRes1=resFiltered3[resFiltered3$log2FoldChange>=1,]
downRes1=resFiltered3[resFiltered3$log2FoldChange<=-1,]

#To view the genes which were upregulated or downregulated
rownames(upRes1)
rownames(upRes2)

rownames(downRes1)
rownames(downRes2)
write.table(resFiltered2,file="all_genelist_output.txt",sep="\t",quote=F,row.names=T, col.names=T)
write.table(upRes1,file="upRes_logfc1_pval0.05.txt",sep="\t",quote=F,row.names=T, col.names=T)
write.table(downRes1,file="DownRes_logfc1_pval0.05.txt",sep="\t",quote=F,row.names=T, col.names=T)