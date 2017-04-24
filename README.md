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
### From bioconda
```
conda install -c bioconda targqc
```

### From PyPI
```
pip install targqc
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


# Venn diagrams for BED files
Build a web-page with size-proportional Venn diagrams for an unlimited set of BED files:
```
bed_venn.py *.bed -o res_dir
```
