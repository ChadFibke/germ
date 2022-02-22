![germ_logo](https://github.com/ChadFibke/germ/blob/main/imgs/germ.png)

## Installation:

Collect germ from github and corresponding singularity container

```
$ git clone git@github.com:ChadFibke/germ.git

$ cd ./germ

$ mkdir container

$ cd container

$ singularity pull --arch amd64 library://chadfibke/workflows/germ:latest .
```

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

## Workflow explained

## Output

```
outdir/
|-- read_summary/
|   |-- test_[12]P.fastq.gz - Fastq files after trimming
|   |-- prefix_[12]P_fastqc.zip - Fastqc summary of both fastq files after trimming
|   |-- prefix_[12]_seqkit.tsv - Seqkit summary of both fastq files after trimming (information on read count, length, quality quantiles)
|
|-- assembly/
    |-- prefix.skesa.fa - Skesa denovo assembly (Draft genome in fasta form)
    |-- prefix.fix.fasta - Corrected draft genome
    |-- prefix._assembly_summary.tsv - Number of contigs, total length of genome, min/avg/max contig length, N50, GC% of genome
    |-- prefix_taxonomy.tsv - Taxonomic identity of draft genome
    |-- prefix.mlst.tsv -The file name, scheme, sequence type and allele IDs
    |-- prefix.abricate.tsv - The file name, start/end (genomic coordinates), strand, gene name, coverage and feature identity
    |--busco/
      |-- short_summary.specific.bacteria_odb10.busco.txt - BUSCO summary for genome completeness     
```
