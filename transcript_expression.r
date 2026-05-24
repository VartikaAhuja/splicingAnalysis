#Filtering overall tumor-specific isoforms (upregulated by logFoldChange >=1 and tumor-specific)

count_data = as.matrix(read.csv("transcript_count_matrix.csv",sep=",",row.names="transcript_id"))
col_data=read.table("samples_info.txt",sep="\t",header=T)
#samples_info file has columns - sample name, type, condition. In tumor samples analysis, condition may be tumor or normal. other columns for stage, site, etc. can be added
rownames(col_data)=col_data$Sample

all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))

dds = DESeqDataSetFromMatrix(countData = count_data,
                               + colData = col_data,
                               + design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "normal") 
rld = rlogTransformation(dds)
library("RColorBrewer") 
library("pheatmap")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition))
colnames(mat) = NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file="sample_pca.pdf")
plotPCA(rld, intgroup=c("condition"))
dev.off()
dds = estimateSizeFactors(dds)
sizeFactors(dds)
dds = DESeq(dds)

res<-results(dds)
head(res)
resOrdered <- res[order(rownames(res)),]
write.table(resOrdered,file="deseq2_all_transcripts_results.txt",sep="\t",quote=F,row.names=T, col.names=T)
resFiltered <- subset(resOrdered,resOrdered$pvalue<0.05)
write.table(resFiltered,file="deseq2_all_transcripts_pval0.05.txt",sep="\t",quote=F,row.names=T, col.names=T)
upRes1=resFiltered[resFiltered$log2FoldChange>=1,]

tumor_cols <- col_data$Sample[col_data$condition=="tumor"]
length(tumor_cols)
normal_cols <- col_data$Sample[col_data$condition=="normal"]
length(normal_cols)
normal_cols
dim(count_data)

new_data=count_data[rownames(upRes1),]
dim(new_data)
#giving a count for each isoform specifying number of tumor samples and number of normal samples with expression upregulated by logFC 1 atleast
count_tumor <- rowSums(new_data[,tumor_cols]>0)
count_normal <- rowSums(new_data[,normal_cols]>0)
head(count_tumor)

#depending on the prevalence, choose n1 and n2 examples - 0.8 and 0.2 for atleast 80% tumor samples and atmost 20% normal samples
#or n1 and n2 as 1 and 0 so 100% tumor samples and none normal samples 
n1=0.8
n2=0.2

threshN <- floor(n2*length(normal_cols))
threshT <- floor(n1*length(tumor_cols))
threshT
threshN
rows_keep <- count_tumor>=threshT & count_normal<=threshN
#indicates the number of isoforms which satisfied the condition
sum(rows_keep)

matr=data.frame(count_tumor,count_normal)
ids <- read.table("tid_geneid.txt",sep="\t",header=T,row.name=1)
#file tid_geneid.txt has transcript ids and their corresponding gene ids in 2 columns, with transcript id as column 1
matr=merge(matr,ids,by.x='row.names',by.y='row.names')      #generates a table with transcript ids, gene ids, and number of tumor and normal samples where they are present
write.table(matr,"upregulated_isoforms_all.tsv",sep="\t",row.names=F,quote=F,col.names=T)