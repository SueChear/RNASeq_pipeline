# RNASeq_pipeline

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
Open count.out contents in spreadsheet, remove first row, keep 'geneid' and 'sample' columns, save as countfile.csv

