# RNA-Seq data analysis pipeline: Reprocessing data from a Batten cerebral organoid study
Gomez-Giro, G., Arias-Fuenzalida, J., Jarazo, J. et al. Synapse alterations precede neuronal damage and storage pathology in a human cerebral organoid model of CLN3-juvenile neuronal ceroid lipofuscinosis. acta neuropathol commun 7, 222 (2019). [https://doi.org/10.1186/s40478-019-0871-7](https://doi.org/10.1186/s40478-019-0871-7)

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

To run multiqc, cd to directory with fastqc outputs.

```
multiqc *fastqc*

```

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
```
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```
Download the annotation gtf file, return to Ensemble Human homepage, click "Download GTF or GFF3". Right click copy link for "homo_sapiens.GRCh38.108.gtf.gz"
At the terminal, type
``` 
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
```
Unzip both Human Reference genome and gtf file.

In STAR-2.7.10b folder, make a new directory named "ref" to house the genome indice file. Ensure that gtf and reference genome file are in STAR folder (not within any folder in STAR)
```
STAR –runMode genomeGenerate –genomeDir ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa –sjdbGTFfile Homo_sapiens.GRCh38.108.gtf –runThreadN 32
```

Create a new folder named 'fastq' in STAR, place all fastq files in 'fastq' folder. Create a 'mapped' folder in STAR.
```
cd fastq
For file in *.fq;do STAR --runMode alignReads --genomeDir ../ref/ --outSAMtype BAM SortedByCoordinate --readFilesIn ${file} --runThreadN 12 --outFileNamePrefix ../mapped/${file};done
```

For MacOS

```
cd fastq
For file in *.fq;do STAR --runMode alignReads --genomeDir ../ref/ --outSAMtype BAM Unsorted --readFilesIn ${file} --runThreadN 12 --outFileNamePrefix ../mapped/${file};done
```

Examine the percentage of mapped reads.

```
cd mapped
cat SRR10116296.fastqLog.final.out
```

Create a directory for bam file named 'bams' in 'mapped' folder and move bam files to 'bams'
```
mkdir bams
cd bams
mv *.bam bams/
```

To check mapping, use samtools to sort and index.(for macOS, the mapped sequence is not sorted)
```
samtools sort usorted.bam > sorted.bam
samtools index sorted.bam
```
Then to view mapped file, go to https://igv.org/app/, load the hg38 ref sequence, followed by the sorted.bam and sorted.bam.bai (index file). Zoom in and out of 
the tracks to view the mapping.

Install subread, use featureCounts to count reads
```
sudo apt-get update
sudo apt-get -y install subread

featureCounts -a Homo_sapiens.GRCh38.108.gtf -o count.out -T 8 mapped/bams/*.bam
```

For MacOS
featureCounts will sort bam files if there were unsorted.

cd to subread folder>bin

```
./featureCounts -a Homo_sapiens.GRCh38.108.gtf -o count.out -T 8 bams/*.bam
```
Open count.out contents in spreadsheet, remove first row, keep 'geneid' and 'sample' columns, save as 'count.csv'

Use multiqc to examine featureCounts summary
```
multiqc count.out.summary
```

## Differential expression
Run all of the following codes in R.

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
Create column called 'condition' which consists of either 'C' for 'control' or 'S' for 'mutant'
```
condition<-factor(c("C","C","C","C","C","C","S","S","S","S","S","S"))

coldata<-data.frame(row.names=colnames(counts),condition)
```
A design formula tells the statistical software the known sources of variation to control for, as well as the factor of interest to tets for during differential expression testing. For example, if you know that batch is a significant source of variation in your data, then batch should be included in your model. The design formula should have all of the factors in your metadata that account for major sources of variation in your data. The last factor entered in the formula should be the condition of interest. For example, here we want to examine the expression differences between the groups, control and sample, then your deisng formula would be design=~condition. The tilde ~ should always proceed your factors and tells DESeq2 to model the counts using the formula. Note the factors included in the design formula need to match the column names in the metadata. If there is another column 'batch' as source of variation, then model formula will be design<-~batch+condition. (https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html) The object class used by DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis is the DESeqDataSet, which will usually be represented in the code here as an object dds.

```
dds<-DESeqDataSetFromMatrix(countData = counts,colData = coldata,design=~condition)
```
Use relevel to specify the reference group:
```
dds$condition<-relevel(dds$condition, ref="C")
```

The standard differential expression analysis steps are wrapped into a single functionm, DESeq. This workflow includes estimating size factors, estimating gene-wise dispersion, fitting curve to gene-wise dispersion estimates, shrinking gene-wise dispersion estimates and glm fit for each gene.
DESeq2 calculates the size factors for each sample using the median of ratios method to normalize the count data. Dispersion is a measure of spread or variability in the data. DESeq2 uses a specific measure of dispersion related to the mean and variance of the data. The DESeq2 dispersion estimates are inversely related to the mean and directly related to variance. Based on this relationship, the dispersion is higher for small mean counts and lower for large mean counts. The dispersion estimates for genes with the same mean will differ only based on their variance. Therefore, the dispersion estimates reflect the variance in gene expression for a given mean value.

```
dds<-DESeq(dds)
```

The variance stablizing transformation VST function calculates a VST from the fitted dispersion-mean relation and then transforms the count data, yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values).
```
vsdata<-vst(dds,blind=F)
```
Plot principal component analysis
```
plotPCA(vsdata,intgroup="condition")
```
<img width="639" alt="pcaplot" src="https://user-images.githubusercontent.com/117556524/206892635-2b65adf1-981c-4d47-9398-eec2ef986c98.PNG">

Plot the per-gene dispersion estimates together with the fitted mean-dispersion relationship. The curve is displayed as a red line in the figure below, which plots the estimate for the expected dispersion value for genes of a given expression strenght. Each black dot is a gene with an associated mean expression level and
maximum likelihood estimation of the dispersion. We expect the dispersion to decrease as the mean of normalized counts increases.

```
plotDispEsts(dds)
```
<img width="941" alt="dispersionestimates" src="https://user-images.githubusercontent.com/117556524/206892888-a1f1a272-0058-40f4-8249-b4e97f86b51d.PNG">


Contrast() to compare between group 'C' and 'S'
```
res<-results(dds, contrast=c("condition","S","C"))
```
Omit NA from data
```
sigs<-na.omit(res)
```
Filter adjusted p value<0.05
```
sigs<-sigs[sigs$padj<0.05,]
```

## Convert ENSEMBL ID to symbol 
Install human database from bioconductor
```
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
```
Convert sigs.df into data frame
```
sigs.df<-as.data.frame(sigs)
```

Use columns to discover which sorts of annotations can be extracted from it
```
columns(org.Hs.eg.db)
```

## Heatmap
Install ComplexHeatmap
```
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
```
Filter basemean>150 and log2fold change>2.5
```
sigs.df<-sigs.df[(sigs.df$baseMean>150) & (abs(sigs.df$log2FoldChange)>2.5),]
```

keytypes are same as columns, keys are values under each keytype. By simply using appropriate argument values with select we can specify what 
keys we want to look up values for (keys), what we want returned back (columns) and the type of keys that we are passing in (keytype), keytype has similar values 
between two datasets, in this case it is "ENSEMBL" where keys are similar to rownames of sigs.df
```
sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype = "ENSEMBL",column="SYMBOL")


mat<-counts(dds, normalized=T)[rownames(sigs.df),]

mat.z<-t(apply(mat,1,scale))

colnames(mat.z)<-rownames(coldata)

Heatmap(mat.z,cluster_rows=T, cluster_columns=T, column_labels=colnames(mat.z), 
        name="Z-score", row_labels=sigs.df[rownames(mat.z),]$symbol, 
        row_names_gp =gpar(fontsize=5, fontfamily="sans", fontface="bold"))
```

<img width="948" alt="heatmap" src="https://user-images.githubusercontent.com/117556524/206892917-dd746a85-7f05-49fe-bb08-02f5ac21d198.PNG">

## Gene ontology enrichment with ClusterProfiler
Install ClusterProfiler
```
BiocManager::install("clusterProfiler")

library(clusterProfiler)
```
Retrieve just the gene ID ENSG (rownames) for log2FC>1 
```
genes_to_test <- rownames(sigs[sigs$log2FoldChange >1,])

```
Retrieve gene ontology terms associated with upregulated genes in 'mutant' samples
```
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)
```
<img width="413" alt="GOtable" src="https://user-images.githubusercontent.com/117556524/206893108-263cae40-9179-4af6-a4df-4d90ff297bd0.PNG">

```
plot(barplot(GO_results, showCategory = 15))
```
<img width="949" alt="GOplot" src="https://user-images.githubusercontent.com/117556524/206893149-9354a7ee-94dd-4d5f-aef4-b3f43f5687e2.PNG">

Retrieve gene ontology terms associated with downregulated genes in 'mutant' samples
```
genes_to_test <- rownames(sigs[sigs$log2FoldChange < -1,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)

plot(barplot(GO_results, showCategory = 15))
```
<img width="949" alt="GOplotdown" src="https://user-images.githubusercontent.com/117556524/206900258-3acf90b3-15ba-4702-b831-f8a8f8d7e2d2.PNG">


## Volcano plot
Install EnhancedVolcano and other packages
```
BiocManager::install('EnhancedVolcano')

install.packages("textshaping")

install.packages("ggrepel")
library(EnhancedVolcano)

```

Using pcutoff=1e-4 and Fold change cutoff=1

```
sigs<-na.omit(res)
sigs<-sigs[sigs$padj<0.05,]
sigs.df<-as.data.frame(sigs)
sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype="ENSEMBL",column="SYMBOL")

EnhancedVolcano(sigs.df, x="log2FoldChange", y="padj", lab=sigs.df$symbol, pCutoff=1e-4, FCcutoff=1)
```
<img width="948" alt="volcano" src="https://user-images.githubusercontent.com/117556524/206893171-a6f3feba-9b6f-40a3-bd03-a41a527ab833.PNG">

To label a few selected genes only

```
selected=c("COL18A1","HOXB9")
EnhancedVolcano(sigs.df, x="log2FoldChange", y="padj", lab=sigs.df$symbol, pCutoff=1e-4, FCcutoff=1, selectLab=selected)
```
<img width="950" alt="volcanolabel" src="https://user-images.githubusercontent.com/117556524/206893192-9860c611-424f-4461-bba7-4d687cdc4080.PNG">

## GSEA analysis 
(https://stephenturner.github.io/deseq-to-fgsea/))

Load required packages.
```
library(fgsea)
library(magrittr)
library(tidyverse)
```

Retrieve gene symbol and test statistic which are what we're intrested in.
```
res_fg2<-sig.df%>%dplyr::select(symbol,stat)%>%na.omit()%>%distinct()%>%
         group_by(symbol)%>%summarize(stat=mean(stat))

```
The fgsea function requires a list of gene sets to check and a named vector of gene-level statistics, where the names should be the same as the gene names in the pa thway list. First let's create our named vector of test statistics. Deframe converts two column data frames to a named vector or list using the first column as name and the second column as value.
```
ranks<-deframe(res_fg2)
head(ranks,20)
```
The gmtpathways() function will take a gmt file you downloaded from MSigDB and turn it into a list. Each element in the list is a character vector of genes in 
the pathway. Load the pathways into a named list.
```
pathways.hallmark<-gmtPathways("h.all.v2022.1.Hs.symbols.gmt")
pathways.hallmark
```

Show the first few pathways, and within those, show only the first few genes
```
pathways.hallmark%>% head() %>% lapply(head)
```

Now run the fgsea algorithm with 1000 permutations
```
fgseaRes<-fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
```
Tidy the results
```
fgseaResTidy<-fgseaRes%>%as_tibble() %>%arrange(desc(NES))
```
Show results in a table
```
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj)%>%DT::datatable()
```

<img width="956" alt="gseatable" src="https://user-images.githubusercontent.com/117556524/206893208-741cc725-753e-4051-902f-b32cade69d23.PNG">
Plot the normalized enrichment scores. Color the bar indicating wheteher or not the pathway was significant.

```
ggplot(fgseaResTidy, aes(reorder(pathway,NES),NES))+ 
  geom_col(aes(fill=padj<0.05))+coord_flip()+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")+
  theme_minimal()
```

<img width="948" alt="hallmarkplot" src="https://user-images.githubusercontent.com/117556524/206893230-41d73776-23fa-4701-9a91-2fdd1c24ab03.PNG">
Show enrichment plot for selected pathway.

```
plotEnrichment(pathway=pathways.hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], ranks)
```

<img width="957" alt="gseaplot" src="https://user-images.githubusercontent.com/117556524/206893295-4cfe326c-6999-4a8a-afa8-8f171a0e31a6.PNG">

PlotGseaTable allows us to plot a summary figure showing multiple pathways.

```
plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj<0.05]],ranks,
              fgseaRes, gseaParam = 0.5)
```

<img width="752" alt="gseasummary" src="https://user-images.githubusercontent.com/117556524/206893308-9229a340-63c9-4d48-8354-5f9c90b889d9.PNG">




