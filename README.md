# germ

## Installation:
`git clone git@github.com:ChadFibke/germ.git`

## Usage:
First add germ to your `$PATH`

Once in your `$PATH`, pass your
`
germ \
    prefix=test \
    r1=/Path/to/test_1.fastq.gz \
    r2=/Path/to/test_2.fastq.gz \
    germ=/Path/to/germ \
    outdir=/Path/to/outdir
`

| Argument | Description |
| ----------- | ----------- |
| prefix | Name of sample. This will be appended to all output |
| r1 | Absolute pathway to fastq file storing read 1 sequences |
| r2 | Absolute pathway to fastq file storing read 2 sequences |
| germ | Absolute pathway to the git cloned germ directory |
| outdir | A directory you would like the results to be written to |
