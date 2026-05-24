# Perform initial quality control with FastQC on raw FASTQ files
# Trim Illumina adapters and low-quality bases using TrimGalore or Cutadapt

## Align samples to the human reference genome
# Build HISAT2 index from hg38.fa (prefix 'hg38' will be used to reference index files)
hisat2-build hg38.fa hg38   #building an index of the reference genome hg38.fa and hg38 is the prefix in order to refer to the index files

# Align paired-end reads; --dta parameter is important, to use hisat aligned bam files for stringtie later
hisat2 -x hg38 --dta -p 10 -1 <sample_R1_val_1>.fq.gz -2 <sample_R2_val_2>.fq.gz -S <sample>.sam

# Convert SAM to BAM, sort, and index for efficient downstream processing
samtools view -bS <sample>.sam -o <sample>.bam -@ 10
samtools sort <sample>.bam -o <sample_sorted>.bam -@ 10
samtools index <sample_sorted>.bam -@ 10

# Assemble transcripts and estimate abundance using StringTie with GENCODE annotation
stringtie -G gencode.v45.basic.annotation.gtf -e -o ./<sample>/<sample>_stringtie.gtf -v -A <sample>_gene_abundance.out -B -p 10 <sample_sorted>.bam


# Generate gene-level counts using bedtools multicov
bedtools multicov -bams <sorted_sample1>.bam <sorted_sample2>.bam ... -bed hg38_Genes.bed > gene_counts.txt
#The file generated doesnot have header information, Add header information to gene count file
sed -i '1s/^/chr\tstart\tend\tgene\t<sample1>\t<sample2>\n/g' gene_counts.txt 

#to generate transcript expression counts:
#create a file containing path to each sample's stringtie output
#for example: samples.csv
#sample1,/path/to/sample1_stringtie.gtf
#sample2,path/to/sample2_stringtie.gtf

# Prepare transcript-level counts using prepDE.py with a sample list
#prepDE.py is a python script utility provided by Stringtie
python prepDE.py -i samples.csv
