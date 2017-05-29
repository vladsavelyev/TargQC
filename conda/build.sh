conda build targqc -c vladsaveliev -c bcbio -c bioconda -c conda-forge --py 2.7
conda build targqc -c vladsaveliev -c bcbio -c bioconda -c conda-forge --py 3.6
cd /Users/vlad/miniconda3/conda-bld
conda convert -p linux-32 osx-64/targqc-1.5.0-py36_2.tar.bz2
conda convert -p linux-64 osx-64/targqc-1.5.0-py36_2.tar.bz2
conda convert -p linux-32 osx-64/targqc-1.5.0-py27_2.tar.bz2
conda convert -p linux-64 osx-64/targqc-1.5.0-py27_2.tar.bz2
anaconda upload osx-64/targqc-1.5.0-py36_2.tar.bz2
anaconda upload osx-64/targqc-1.5.0-py27_2.tar.bz2
anaconda upload linux-32/targqc-1.5.0-py36_2.tar.bz2
anaconda upload linux-32/targqc-1.5.0-py27_2.tar.bz2
anaconda upload linux-64/targqc-1.5.0-py36_2.tar.bz2
anaconda upload linux-64/targqc-1.5.0-py27_2.tar.bz2

