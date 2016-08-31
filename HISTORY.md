# 1.2.0
- Downsampling to 5% by default instead of 500k read pairs

# 1.1.0
- Moved code from AZ Reporting Suite.
- Decoupled from SGE using iPython-cluster-helper (https://github.com/roryk/ipython-cluster-helper).
- New BED annotation, using Ensembl
- New region coverage reports
- Support FastQ processing with downsampling
- Add nosetests
- Add Travis CI support
- Dockerized
- Added MultiQC module
- Much less dependencies, easier to install
- Setup.py and uploaded PyPI (available via pip)