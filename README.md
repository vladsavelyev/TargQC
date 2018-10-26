# TargQC - target capture coverage QC tool

[![Anaconda-Server Badge](https://anaconda.org/vladsaveliev/targqc/badges/installer/conda.svg)](https://conda.anaconda.org/vladsaveliev)
[![Build Status](https://travis-ci.org/vladsaveliev/TargQC.svg?branch=master)](https://travis-ci.org/vladsaveliev/TargQC)

## Input

- BAM file(s) (or FastQ files).
- BED file (optional).

## Output

- `summary.html` – sample-level coverage statistics and plots.
- `summary.tsv` – sample-level coverage, parsable version
- `regions.tsv` – region-level coverage statistics.

To view `regions.tsv` contents in a nice aligned table, you can use `tsv` command, e.g.:

```
tsv regions.tsv
chr    start     end       size  gene          exon        strand  feature  biotype                             transcript       trx_overlap  exome_overlap  cds_overlap  avg_depth  at1x     at5x     at10x     at20x    at50x     at100x   at250x   at500x   at1000x  at5000x  at10000x  at50000x
chr21  9907173   9907501   328   TEKT4P2       3           -       capture  transcribed_unprocessed_pseudogene  ENST00000400754  95.1%        95.1%          0%           115.137    100      100      100       100      100       66.7683  0        0        0        0        0         0
chr21  10863011  10863111  100   IGHV1OR21-1   2           +       capture  IG_V_gene                           ENST00000559480  56.0%        56.0%          53.0%        106.93     100      100      100       100      100       100      0        0        0        0        0         0
chr21  10910285  10910401  116   TPTE          22          -       capture  protein_coding                      ENST00000361285  100.0%       80.2%          80.2%        36.3621    100      100      100       100      0         0        0        0        0        0        0         0
```

Or you can use bioawk (`conda install -c bioconda bioawk`) to query certain columns with:

```
cat ./tests/gold_standard/bed3/syn3-tumor/regions.tsv | bioawk -tc hdr '{ print $chr, $start, $end, $gene, $avg_depth, $at100x }' | tsv
chr    start     end       gene        avg_depth  at100x
chr21  9907173   9907501   TEKT4P2     115.137    66.7683
chr21  10863011  10863111  IGHV1OR21-1 106.93     100    
chr21  10910285  10910401  TPTE        36.3621    0      
```

## Columns explanation
 
- `trx_overlap`, `exome_overlap`, `cds_overlap` - percentage of the region that overlaps transcripts, exons, or CDS (coding regions) in Ensembl database, correspondingly.

- `at1x`, `at5x`, etc - percentage of the region that is covered at least at this depth (1x, 5x) by reads, excluding duplicated, quality-failed, secondary and supplemetrary aligned reads.


## Installation

### From conda

Install miniconda if not already:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh   # linux
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh  # macos
bash miniconda.sh -b -p ./miniconda
unset PYTHONPATH
export PATH=$(pwd)/miniconda/bin:$PATH
conda activate base  # or create an environment with `conda env create -n targqc`
```

Install targqc:

```
conda install -c vladsaveliev -c bioconda -c conda-forge targqc
```


### From source

This approach assumes you have dependencies like qualimap, bedtools, sambamba already installed on your system and available in PATH.

```
git clone --recursive https://github.com/vladsaveliev/TargQC
cd TargQC
virtualenv venv_targqc && source venv_targqc/bin/activate  # optional, but recommended if you are not an admin
pip install --upgrade pip setuptools
pip install -r requirements.txt
python setup.py install
```

## Usage

```
targqc *.bam --bed target.bed -g hg19 -o targqc_results
```

The results will be written to `targqc_results` folder.

The BED file may be omitted. In this case statistics reported will be based of off the whole genome.

The accepted values for `-g` are `hg19`, `hg38`, or a full path to any indexed reference fasta file:

```
targqc *.bam --bed target.bed -g /path/to/genomes/some_genome.fa -o targqc_results
```

When running from BAMs, only the `.fai` index is used, and the fasta file itself can be non-existent.

Instead of the BAM files, input FastQ are also allowed. The reads will be aligned by BWA to the reference 
genome specified by `--bwa-prefix` (unless `-g` is already a fasta path bwa-indexed).

```
targqc *.fastq --bed target.bed -g hg19 -o targqc_results --bwa-prefix /path/to/ref.bwa
```

Option `--downsample-to <N>` (default value `5e5`) specifies the number of 
read pairs will be randomly selected from each input set. This feature allows to quickly estimate approximate 
coverage quality before full alignment. To turn downsampling off and align all reads, set `--downsample-to off`.


## Parallel running

### Threads

Run using 3 threads:

```
targqc *.bam --bed target.bed -g hg19 -o targqc_results -t 3
```

### Cluster

Run using 3 jobs, using SGE scheduler, and queue "queue":

```
targqc *.bam --bed target.bed -g hg19 -o targqc_results -t 3 -s sge -q batch.q -r pename=smp
```

If the number of samples is higher than the requested number of jobs, the processes within job will be additionally parallelized using threads, so the full number of occupied cores will equal the number of requested threads (-t)

Other supported schedulers: Platform LSF ("lsf"), Sun Grid Engine ("sge"), Torque ("torque"), SLURM ("slurm") (see details at https://github.com/roryk/ipython-cluster-helper)


# BED Annotation

A tool that assings gene names to regions in a BED file based on Ensembl genomic features overlap.

[![Anaconda-Server Badge](https://anaconda.org/vladsaveliev/bed_annotation/badges/installer/conda.svg)](https://conda.anaconda.org/vladsaveliev)

UPD: Moved into a separate repository [https://github.com/vladsaveliev/bed_annotation](https://github.com/vladsaveliev/bed_annotation)


# Venn diagrams for BED files

Build a web-page with size-proportional Venn diagrams for an unlimited set of BED files:

UPD: Moved into a separate repository [https://github.com/vladsaveliev/Venn](https://github.com/vladsaveliev/Venn)

