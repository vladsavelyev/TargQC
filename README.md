[![Anaconda-Server Badge](https://anaconda.org/vladsaveliev/targqc/badges/installer/conda.svg)](https://conda.anaconda.org/vladsaveliev)
[![Build Status](https://travis-ci.org/vladsaveliev/TargQC.svg?branch=master)](https://travis-ci.org/vladsaveliev/TargQC)

# TargQC - target capture coverage QC tool

## Input

- BAM file(s) (or FastQ files).
- BED file (optional).

## Output

- `summary.html` – sample-level coverage statistics and plots.
- `summary.tsv` – sample-level coverage, parsable version
- `regions.txt` – region-level coverage statistics.

## Installation

### From conda

```
conda install -c vladsaveliev targqc
```

### From source

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


# BED_Annotation - a tool for BED file annotation with Ensemble gene names

[![Anaconda-Server Badge](https://anaconda.org/vladsaveliev/bed_annotation/badges/installer/conda.svg)](https://conda.anaconda.org/vladsaveliev)

The package provides a script named `annotate_bed.py` that annotates a BED file regions with gene symbols, based on Ensembl data.

### Usage

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed
``` 

Script checks each region against the Ensembl genomic features database, and writes a BED file in a standardized format with a gene symbol, strand and exon rank in 4-6th columns:

`INPUT.bed`:

```
chr1    69090   70008
chr1    367658  368597
```

`OUTPUT.bed`:

```
chr1    69090   70008   OR4F5   1       +
chr1    367658  368597  OR4F29  1       +
```

#### Transcripts order

The piority for choosing transcripts for annotation is the following:
- Overlap % with transcript
- Overlap % with CDS
- Overlap % with exons
- Biotype (`protein_coding` > others > `*RNA` > `*_decay` > `sense_*` > `antisense` > `translated_*` > `transcribed_*`)
- TSL (1 > NA > others > 2 > 3 > 4 > 5)
- Presence of a HUGO gene symbol
- Is cancer canonical
- Transcript size

#### Extended annotation

Use `--extended` option to report extra columns with details on features, biotype, overlapping transcripts and overlap sizes:

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended
```

`OUTPUT.bed`:

```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon  Strand  Feature Biotype         Ensembl_ID      TSL HUGO    Tx_overlap_% Exon_overlaps_% CDS_overlaps_% Ori_Fields
chr1    69090   70008   OR4F5   1     +       capture protein_coding  ENST00000335137 NA  OR4F5   100.0        100.0           99.7
chr1    367658  368597  OR4F29  1     +       capture protein_coding  ENST00000426406 NA  OR4F29  100.0        100.0           99.7
```

#### Ambuguous annotations

Regions may overlap mltiple genes. The `--ambiguities` controls how the script resolves such ambiguities

- `--ambiguities all` -- report all reliable overlaps (in order in the "priority" section, see above)
- `--ambiguities all_ask` -- stop execution and ask user which annotation to pick
- `--ambiguities best_all` (default) -- find the best overlap, and if there are several equally good, report all (in terms of the "priority" above)
- `--ambiguities best_ask` -- find the best overlap, and if there are several equally good, ask user
- `--ambiguities best_one` -- find the best overlap, and if there are several equally good, report any of them

Note that the first 4 options might output multiple lines per region, e.g.:

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended --ambiguities best_all
```

`OUTPUT.bed`:

```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000412321 NA      OR4F29  100.0   100.0   100.0
```

#### Other options

- `--coding-only`: take only the features of type `protein_coding` for annotation
- `--high-confidence`: annotate with only high confidence regions (TSL is 1 or NA, with HUGO symbol, total overlap size > 50%)
- `--canonical`: use only canonical transcripts to annotate (which to the most part means the longest transcript, by SnpEff definition)
- `--short`: add only the 4th "Gene" column (outputa 4-col BED file instead of 6-col)
- `--output-features`: good for debugging. Under each BED file region, also output Ensemble featues that were used to annotate it
 

# Venn diagrams for BED files

Build a web-page with size-proportional Venn diagrams for an unlimited set of BED files:

```
bed_venn.py *.bed -o res_dir
```
