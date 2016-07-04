# Alignment target coverage analysis tool 

## Input
- BAM file(s) (or FastQ files).
- BED file (optional).

## Output
- `summary.html` – sample-level coverage statistics and plots.
- `summary.tsv` – sample-level coverage, parsable version
- `regions.txt` – region-level coverage statistics.

## Installation
```
git clone --recursive https://github.com/vladsaveliev/TargQC
cd TargQC
python setup.py install
```

## Usage
```
targqc *.bam --bed target.bed -g hg19 -o targqc_results
```
The results will be written to `targqc_results` folder.

The BED file may be omitted. In this case statistics reported will be based of off the whole genome.

## FastQ and downsampled coverage
Instead of the BAM files, input FastQ are also allowed. The reads will be aligned by BWA to the reference genome specified by `--bwa-prefix`. Option `--downsample-to N` (default value `5e5`) specifies the number ofread pairs will be randomly selected from each input set. This feature allows to quickly estimate approximate coverage quality before full alignment. To turn downsampling off and align all reads, set `--downsample-to off`.
```
targqc *.fastq --bed target.bed -g hg19 -o targqc_results --bwa-prefix /path/to/ref.bwa
```

## Parallel running
### Threads
Run using 3 threads:
```
targqc *.bam --bed target.bed -g hg19 -o targqc_results -t 3
```
### Cluster
Run using 3 jobs, using SGE scheduler, and queue "queue":
```
targqc *.bam --bed target.bed -g hg19 -o targqc_results -t 3 -s sge -q queue -r pename=smp
```
If the number of samples is higher than the requested number of jobs, the processes within job will be additionally parallelized using threads, so the full number of occupied cores will equal the number of requested threads (-t)

Other supported schedulers: Platform LSF ("lsf"), Sun Grid Engine ("sge"), Torque ("torque"), SLURM ("slurm") (see details at https://github.com/roryk/ipython-cluster-helper)
