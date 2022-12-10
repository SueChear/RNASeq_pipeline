# RNA-Seq data analysis pipeline tutorial with an example of repurposing data from a Batten cerebral organoid study

GEO datasets are downloaded using SRA toolkit. To download SRA toolkit, follow instructions at https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit. Download tar file for Ubuntu from NCBI.
```
$ wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```
Extract contents of the tar file
```
tar -vxzf sratoolkit.tar.gz
```
Append the path to the binaries to your PATH environment variable. Path environment variable specifies the directories to be searched to find a command, here the path environment variable is 'bin' which usually contains all executable files.
If sratoolkit.3.0.0-ubuntu64 folder is at home/user/rnaseq,

```
cd rnaseq
export PATH=$PATH:$PWD/sratoolkit.3.0.0-ubuntu64/bin
```
Verify that the binaries will be found by the shell:
```
which fastq-dump
```
This should produce output similar to
```
/home/user/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump
```
Define download path and other tool settings with SRA-toolkit configuration
```
vdb-config --interactive
```
This brings you to a window, use "Tab" on keyboard to select "Main". Ensure "Enable remote access" is selected, then go to "Cache", ensure "Enable local file caching" is selected, then go to "Location of user repository", select folder to download SRA files to. Save changes and exit.

To download SRR file
```
fasterq-dump --split-files SRR********
```
To show first 20 lines of fastq file
```
cat SRR********.fastq | head -n 20
```

## Quality control of fastq files
Install fastqc using apt-get, update apt database with apt-get using the following command.
```
sudo apt-get update
```
Install fastqc
```
sudo apt-get -y install fastqc
```
Run fastqc on fastq files. Move all fastq files for qc into one folder, eg to folder named "fastqc"
```
cd fastqc
fastqc *.fastq
```
For differential expressed gene analysis, trimming of low quality sequences and adapter is not necessary.

## Generating genome indices
Download STAR
```
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b
sudo apt install rna-star
STAR --version
```
To download human genome, go to Ensemble Human, click "Download DNA sequence (FASTA)". Right click copy link for "homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
At the terminal, type
```wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```
Download the annotation gtf file, return to Ensemble Human homepage, click "Download GTF or GFF3". Right click copy link for "homo_sapiens.GRCh38.108.gtf.gz"
At the terminal, type
``` wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
```
Unzip both Human Reference genome and gtf file.

In STAR-2.7.10b folder, make a new directory named "ref" to house the genome indice file. Ensure that gtf and reference genome file are in STAR folder (not within any folder in STAR)
```
STAR –runMode genomeGenerate –genomeDir ref/ --genomeFastaFiles Homo_sapies.GRCh38.dna_sm.primary_assembly.fa –sjdbGTFfile Homo_sapiens.GRCh38.108.gtf –runThreadN 32
```

Create a new folder named 'fastq' in STAR, place all fastq files in 'fastq' folder. Create a 'mapped' folder in STAR.
```
cd fastq
For file in *.fq;do STAR --runMode alignReads --genomeDir ../ref/ --outSAMtype BAM SortedByCoordinate --readFilesIn ${file} --runThreadN 12 --outFileNamePrefix ../mapped/${file};done
```

Create a directory for bam file named 'bams' in 'mapped' folder and move bam files to 'bams'
```
mkdir bams
cd bams
mv *.bam bams/
```
Install subread, use featureCounts to count reads
```
sudo apt-get update
sudo apt-get -y install subread

featureCounts -a Homo_sapiens.GRCh38.108.gtf -o count.out -T 8 mapped/bams/*.bam
```
Open count.out contents in spreadsheet, remove first row, keep 'geneid' and 'sample' columns, save as 'count.csv'

##Differential expression
All following code should be run in R.

Install the required packages:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
```

Import count data
```
counts<-read.delim("count.csv",header=T, row.names=1, sep=",")
head(counts)
```

Filter counts >50
```
counts<-counts[which(rowSums(counts)>50),]
```
Create column called 'condition' which consists of either 'control' or 'mutant'
```
condition<-factor(c("C","C","C","C","C","C","S","S","S","S","S","S"))

coldata<-data.frame(row.names=colnames(counts),condition)
```
The design, as shown below indicates how to model the samples, if we want to measure the effect of the condition, then
design=~condition where condition consists of columns in coldata.

```
dds<-DESeqDataSetFromMatrix(countData = counts,colData = coldata,design=~condition)
```
The standard differential expression analysis steps are wrapped into a single functionm, DESeq.
```
dds<-DESeq(dds)
```

The variance stablizing transformation VST function calculates a VST from the fitted dispersion-mean relation and then transforms the count data, yielding a matrix of values which are now approximately hommoskedastic (having constant variance along the range of mean values).
```
vsdata<-vst(dds,blind=F)
```
Plot principal component analysis
```
plotPCA(vsdata,intgroup="condition")
```
Plot dispersion estimates
```
plotDispEsts(dds)
```

```
res<-results(dds, contrast=c("condition","S","C"))

sigs<-na.omit(res)

sigs<-sigs[sigs$padj<0.05,]
```

##Convert ENSEMBL ID to symbol 
Install human database from bioconductor
```
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
```

```
sigs.df<-as.data.frame(sigs)
```

Use columns to discover which sorts of annotations can be extracted from it
```
columns(org.Hs.eg.db)
```

keytypes same as columns, keys are values under each keytype. By simply using appropriate argument values with select we can specify what 
keys we want to look up values for (keys), what we want returned back (columns) and the type of keys that we are passing in (keytype), keytype has similar values 
between two datasets, in this case it is "ENSEMBL" where keys are similar to rownames of sigs.df

```
keytypes(org.Hs.eg.db)
sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype="ENSEMBL", column="SYMBOL")
```
##Heatmap
Install ComplexHeatmap
```
BiocManager::install("ComplexHeatmap")
```
Filter basemean and log2fold change
```
sigs.df<-sigs.df[(sigs.df$baseMean>150) & (abs(sigs.df$log2FoldChange)>2.5),]
sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype = "ENSEMBL",column="SYMBOL")


mat<-counts(dds, normalized=T)[rownames(sigs.df),]

mat.z<-t(apply(mat,1,scale))

coldata

colnames(mat.z)<-rownames(coldata)
mat.z


Heatmap(mat.z,cluster_rows=T, cluster_columns=T, column_labels=colnames(mat.z), 
        name="Z-score", row_labels=sigs.df[rownames(mat.z),]$symbol)
```

##ClusterProfiler
Install ClusterProfiler
```
BiocManager::install("clusterProfiler")

library(clusterProfiler)
```
Retrieve just the gene ID ENSG (rownames)... for log2FC>0.5
```
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5,])

```
Retrieve gene ontology terms
```
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)

plot(barplot(GO_results, showCategory = 15))
```
##Volcano plot
Install EnhancedVolcano and other packages
```
BiocManager::install('EnhancedVolcano')

install.packages("textshaping")

install.packages("ggrepel")
```
Using pcutoff=1e-4
```
EnhancedVolcano(sigs.df, x="log2FoldChange", y="padj", lab=sigs.df$symbol, pCutoff=1e-4, FCcutoff=1)
```
To label a few selected genes only
```
selected=c("COL18A1","HOXB9")
EnhancedVolcano(sigs.df, x="log2FoldChange", y="padj", lab=sigs.df$symbol, pCutoff=1e-4, FCcutoff=1, selectLab=selected)
```
##GSEA analysis
```
res<-na.omit(res)
res<-res[res$baseMean>50,]

res1<-res[order(-res$stat),]
res1

gene_list<-res1$stat
gene_list
names(gene_list)<-rownames(res1)
gene_list

gse<-gseGO(gene_list, ont="BP", keyType = "ENSEMBL", 
           OrgDb = "org.Hs.eg.db", eps=1e-300)

as.data.frame(gse)

gseaplot(gse, geneSetID = 50)

```




