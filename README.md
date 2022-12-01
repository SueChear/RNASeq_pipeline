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

##Quality control of fastq files
Install fastqc using apt-get, update apt database with apt-get using the following command.
```
sudo apt-get update
```
Install fastqc
```
sudo apt-get -y install fastqc
```
