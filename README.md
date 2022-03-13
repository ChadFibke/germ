<img src="https://github.com/ChadFibke/germ/blob/main/imgs/logo.png" width="300" height="150">

## Introduction 

A containerized nextflow pipeline designed to annotate prokaryote genomes. Germ accepts Illumina paired-end reads, performs read processing, denovo genome assembly and annotation of the assembled draft genome. This workflow is a modified, packaged version of the analyses used in a previous [paper]( https://academic.oup.com/ofid/article/6/11/ofz431/5583892?login=false) interested in characterizing Extraintestinal Pathogenic Escherichia coli genomes in an epidemiological context.


## Installation:

Collect germ from github and corresponding singularity container

```
$ git clone git@github.com:ChadFibke/germ.git

$ cd ./germ

$ mkdir container

$ cd container

$ singularity pull --arch amd64 library://chadfibke/workflows/germ:latest .
```
## Workflow

<img src="https://github.com/ChadFibke/germ/blob/main/imgs/workflow.png" width="900" height="500">

1. Paired-end reads are trimmed using [trimmomatic]( https://github.com/usadellab/Trimmomatic). Trimmomatic is used to remove adapter sequences present from adapter read-through. This step also removes low quality bases (phred score < 5) from the beginning and end of each read, and finally uses a sliding window (10bp) approach to remove the 3’ section of the read once an average base of <= 19 is found. After trimming, reads with lengths >= 45bp are kept for sequential analyses.

2. Both [seqkit]( https://bioinf.shenwei.me/seqkit/) and [fastqc]( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are used to summarize the trimmed reads. The total number of reads, read-length summaries, phred-score summaries and all fastqc basic statistics are collected from each read file.

3. The trimmed reads are then assembled into contigs using [skesa]( https://github.com/ncbi/SKESA). The assembly’s contig bases are then corrected using [pilon]( https://github.com/broadinstitute/pilon). Briefly, this correction process will re-align the trimmed reads to the assembled contigs and only retain high quality alignments (mapq >= 20) and bases (phred >=20). Next, the retained alignment bases are summarized and compared to the assembly’s nucleotide at each genomic site across the assembly. Positions across the assembly will be changed to an alternative nucleotide with greater support than the originally placed nucleotide found in the assembly. The contigs are then removed if they have an average coverage less than 10X and a length < 500bp.

4. The corrected assembly is then summarized and validated. First the corrected assembly’s length, contig number, length summary, N50, GC%, and more are summarized. The assembly’s taxonomic origin is also estimated using [mash]( https://github.com/marbl/Mash) and [GTDB]( https://gtdb.ecogenomic.org/). The final validation is to estimate the completeness of the assembly with [busco]( https://busco.ezlab.org/), which identified the proportion of expected single copy orthologs are found in the assembly.

## Usage:
First add germ to your `$PATH`

Once in your `$PATH`, germ can be used with:
```
germ \
    prefix=test \
    r1=/Path/to/test_1.fastq.gz \
    r2=/Path/to/test_2.fastq.gz \
    germ=/Path/to/germ \
    outdir=/Path/to/outdir
```

| Argument | Description |
| ----------- | ----------- |
| prefix | Name of sample. This will be appended to all output |
| r1 | Absolute pathway to fastq file storing read 1 sequences |
| r2 | Absolute pathway to fastq file storing read 2 sequences |
| germ | Absolute pathway to the git cloned germ directory |
| outdir | A directory you would like the results to be written to |


## Output

```
outdir/
|-- read_summary/
|   |-- prefix_[12]P.fastq.gz - Fastq files after trimming
|   |-- prefix.fastqc[12].summary.txt - Fastqc summary for both fastq files after trimming
|   |-- prefix.seqkit[12].tsv - Seqkit summary of both fastq files after trimming (information on read count, length, quality quantiles)
|
|-- assembly/
    |-- prefix.skesa.fa - Skesa denovo assembly (Draft genome in fasta form)
    |-- prefix.fix.fasta - Corrected draft genome
    |-- prefix.assembly.summary.tsv - Number of contigs, total length of genome, min/avg/max contig length, N50, GC% of genome
    |-- prefix.taxonomy.tsv - Taxonomic identity of draft genome
    |-- prefix.mlst.tsv -The file name, scheme, sequence type and allele IDs
    |-- prefix.abricate.tsv - The file name, start/end (genomic coordinates), strand, gene name, coverage and feature identity
    |-- prefix.busco.csv - BUSCO summary for genome completeness     
```
