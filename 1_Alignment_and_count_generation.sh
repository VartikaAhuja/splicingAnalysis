#The initial data cleaning comprises usage of FastQC to check data quality of .fastq files
#This is followed by trimming of adaptors (used for Illumina sequencing) and other low quality regions using tools like trimgalore or cutadapt

## Alignment of samples to reference genome
#Install Hisat2
#Download human reference genome files (hg38.fa as the latest one while writing the code)

hisat2-build hg38.fa hg38   #building an index of the reference genome hg38.fa and hg38 is the prefix in order to refer to the index files
hisat2 -x hg38 --dta -p 10 -1 <sample_R1_val_1>.fq.gz -2 <sample_R2_val_2>.fq.gz -S <sample>.sam
##--dta parameter is important, to use hisat aligned bam files for stringtie later
samtools view -bS <sample>.sam -o <sample>.bam -@ 10
samtools sort <sample>.bam -o <sample_sorted>.bam -@ 10
samtools index <sample_sorted>.bam -@ 10

#Download stringtie and gtf file for gene and transcript information
stringtie -G gencode.v45.basic.annotation.gtf -e -o ./<sample>/<sample>_stringtie.gtf -v -A <sample>_gene_abundance.out -B -p 10 <sample_sorted>.bam


#to generate gene expression counts:
bedtools multicov -bams <sorted_sample1>.bam <sorted_sample2>.bam ... -bed hg38_Genes.bed > gene_counts.txt
#The file generated doesnot have header information, so we add sample headers
sed -i '1s/^/chr\tstart\tend\tgene\t<sample1>\t<sample2>\n/g' gene_counts.txt 

#to generate transcript expression counts:
#create a file containing path to each sample's stringtie output
#for example: samples.csv
#sample1,/path/to/sample1_stringtie.gtf
#sample2,path/to/sample2_stringtie.gtf

python prepDE.py -i samples.csv