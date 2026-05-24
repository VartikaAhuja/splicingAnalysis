#DESeq2 to shortlist upregulated splice factors and filter splice factor-specific (downstream genes) tumor-specific transcripts

library(DESeq2)
getwd()
count_data=read.table("gene_counts.txt",sep="\t",header=TRUE,row.names="gene")
#the first 3 columns for gene count data have chromosome, start and end information for genes
count_data = count_data[,-c(1,2,3)]
head(count_data)
col_data=read.table(file = "samples_info.txt", header = T, sep = "\t")
#samples_info file has columns - sample name, type, condition. In tumor samples analysis, condition may be tumor or normal. other columns for stage, site, etc. can be added
rownames(col_data) = col_data$Sample
all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))

dds = DESeqDataSetFromMatrix(countData = count_data,
                               + colData = col_data,
                               + design = ~ condition)      #features like stage, site and batch if available, they can be used here
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "normal") 
rld = rlogTransformation(dds)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition))
colnames(mat) = NULL
dds = estimateSizeFactors(dds)
sizeFactors(dds)
dds = DESeq(dds)

res<-results(dds)
head(res)
resOrdered <- res[order(rownames(res)),]
resFiltered <- subset(resOrdered,resOrdered$pvalue<0.05)

#To check the upregulated genes among a pre-defined set of splicing factors for study; where genelist has those splicing factors list
dat<- read.table("genelist",sep="\t",header=F,col.names="gene")
resFiltered2 <- subset(resOrdered,rownames(resOrdered) %in% dat$gene)
dim(resFiltered)
dim(resFiltered2)
resFiltered3 <- subset(resFiltered2,resFiltered2$pvalue<0.05)
dim(resFiltered3)

upRes1=resFiltered3[resFiltered3$log2FoldChange>=1,]        #we obtained the upregulated splicing factors from our list of genes of focus

#suppose, upRes1 has a splicing factor - CELF6, so further, we make a file listing the downstream genes spliced by CELF6 and use it here

rm(count_data,dds,keep,rld,distsRL,mat,res,resOrdered,resFiltered,dat,resFiltered2,resFiltered3)

#DESeq2 for transcripts

count_data = as.matrix(read.csv("transcript_count_matrix.csv",sep=",",row.names="transcript_id"))
dim(count_data)
rownames(col_data) = col_data$Sample
all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))

dds = DESeqDataSetFromMatrix(countData = count_data,
                               + colData = col_data,
                               + design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "normal") 
rld = rlogTransformation(dds)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition))
colnames(mat) = NULL
dds = estimateSizeFactors(dds)
sizeFactors(dds)
dds = DESeq(dds)

res<-results(dds)
resOrdered <- res[order(rownames(res)),]
resFiltered <- subset(resOrdered,resOrdered$pvalue<0.05)

#loading the downstream genes list for celf6
sf <- read.table("celf6_genes_literature.txt",sep="\t",header=F,col.names="gene")
#file tid_geneid has information of transcript ids in column 1 and gene ids in column 2
ids <- read.table("tid_geneid.txt",sep="\t",header=T,row.name=1)

sf_2<-subset(ids,ids$gene_name %in% sf)     #filters transcript ids which belong to genes affected by celf6

sfList <- subset(resOrdered,rownames(resOrdered) %in% rownames(sf_2))
dim(sfList)
head(sfList)
#to give gene ids to DESeq2 results for transcripts
for (i in 1:nrow(sfList)) {
sfList$gene[i]=sf_2$gene_name[which(rownames(sf_2)==rownames(sfList)[i])]
}        

#transcripts with logFC>=1 and pvalue<=0.05
sf_up=subset(sfList,((sfList$log2FoldChange>=1)&(sfList$pvalue<0.05)))
dim(sf_up)
sf_down=subset(sfList,((sfList$log2FoldChange<=-1)&(sfList$pvalue<0.05)))
write.table(sfList,"celf6_all_Dtranscripts_expression.txt",sep="\t",quote=F,row.names=T,col.names=T)        #Dtranscripts - downstream transcripts
write.table(sf_up,"celf6_Dtranscripts_upreg_logFC1_transcripts.txt",sep="\t",quote=F,row.names=T,col.names=T)
write.table(sf_down,"celf6_Dtranscripts_downreg_logFC1_transcripts.txt",sep="\t",quote=F,row.names=T,col.names=T)

sum(col_data$condition=="tumor")
sum(col_data$condition=="normal")

#subset of overexpressed transcript counts of celf6 downstream genes from overall transcript count table
counts_filtered <- subset(count_data,rownames(count_data) %in% rownames(sf_up))
dim(sf_up)
dim(counts_filtered)

write.table(counts_filtered,"sf_celf6_upreg_isoforms_logFC1_countmatrix.txt",sep="\t",quote=F,row.names=T,col.names=T)
tumor_cols <- col_data$Sample[col_data$condition=="tumor"]
length(tumor_cols)
tumor_cols
normal_cols <- col_data$Sample[col_data$condition=="normal"]
length(normal_cols)
normal_cols
#generates the prevalence count in tumor and normal samples
count_tumor <- rowSums(counts_filtered[,tumor_cols]>0)
count_normal <- rowSums(counts_filtered[,normal_cols]>0)
n_tumor <- length(tumor_cols)
n_normal <- length(normal_cols)

#atleast n1 samples in tumor and atmost n2 samples in normal
n1=0.8
n2=0.2

threshN <- floor(n2*n_normal)
threshT <- floor(n1*n_tumor)
threshT
threshN

rows_keep <- count_tumor>=threshT & count_normal<=threshN
sum(rows_keep)
samplecounts_isoforms_sfcelf6 <- cbind(count_tumor,count_normal)
#rows_keep list the resultant final isoforms
rows_keep2 = cbind(sfList[names(rows_keep[rows_keep]),],cbind(count_tumor[names(rows_keep[rows_keep])],count_normal[names(rows_keep[rows_keep])]))
write.table(rows_keep2,"celf6_Dtranscripts_80T_20N.tsv",sep="\t",quote=F,row.names=T,col.names=T)
write.table(samplecounts_isoforms_sfcelf6,"IsoformsPrevalence_TumorNormal.tsv",sep="\t",row.names=T,col.names=T,quote=F)