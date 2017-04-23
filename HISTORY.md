# 1.4.4
- Created bioconda package
- Add script to build Venn diagram for a set of BED files:
```
bed_venn.py *.bed -o res_dir
```

# 1.4.1
- Support custom genomes. Provide a path to indexed fasta file with `-g`:
```
targqc *.bam --bed target.bed -o results -g /path/to/some_genome.fa
```
When running from BAMs, only the `.fai` index is used, and the fasta file itself can be non-existent.
When running from fastq, bwa indexes are also required.

# 1.4
- Update BED file annotation package.

# 1.3
- BED annotation: --seq2c mode.

# 1.2
- Downsampling to 5% by default instead of 500k read pairs.

# 1.1
- Moved code from AZ Reporting Suite
- Decoupled from SGE using iPython-cluster-helper (https://github.com/roryk/ipython-cluster-helper)
- New BED annotation, using Ensembl
- New region coverage reports
- Support FastQ processing with downsampling
- Add nosetests
- Add Travis CI support
- Dockerized
- Added MultiQC module
- Much less dependencies, easier to install
- Setup.py and uploaded PyPI (available via pip)