#!/usr/bin/env python
import sys

import os
from os.path import join, isfile, abspath, dirname

import pip
from setuptools import setup, find_packages


print('Upgrading pip and setuptools...')
try:
    pip.main(['install', '--upgrade', 'setuptools', 'pip'])
except StandardError:
    sys.stderr.write('Cannot update pip and setuptools, that might cause errors '
                     'during the following intallation\n')


try:
    from ngs_utils import setup_utils
    from ngs_utils.file_utils import which
    from ngs_utils.setup_utils import run_cmdl
except ImportError:
    print('Installing NGS_Utils...')
    pip.main(['install', 'git+git://github.com/vladsaveliev/NGS_Utils.git'])
    from ngs_utils import setup_utils
    from ngs_utils.file_utils import which
    from ngs_utils.setup_utils import run_cmdl


# if not all(which(tool) for tool in ['bedtools', 'sambamba', 'bwa']):
#     conda_path = join(os.getcwd(), 'anaconda')
#     from sys import platform
#     if platform == "darwin":
#         run_cmdl('wget http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh')
#         run_cmdl('bash Miniconda2-latest-MacOSX-x86_64.sh -b -p ' + conda_path)
#     else:
#         run_cmdl('wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh')
#         run_cmdl('bash Miniconda-latest-Linux-x86_64.sh -b -p ' + conda_path)
#     run_cmdl('PATH=' + join(conda_path, 'bin') + ':$PATH conda install --yes -c bioconda -c conda-forge htslib=1.3 bedtools sambamba bwa -q')


name = 'TargQC'
script_name = 'targqc'
package_name = 'targqc'

version = setup_utils.init(name, package_name, __file__)


setup(
    name=name,
    version=version,
    author='Vlad Saveliev and Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='Genome capture target coverage evaluation tool',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/TargQC',
    download_url='https://github.com/vladsaveliev/TargQC/releases',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        package_name: [
            'bedops/bedops_*',
            'qualimap/*/qualimap',
            'qualimap/*/qualimap.jar',
            'qualimap/*/lib/*.jar',
            'qualimap/*/scripts/*.jar',
            'qualimap/*/species/*.jar',
            'picard/picard/*.jar',
            'gender/*.bed',
            'bwa/bwa',
        ],
        'ensembl': [
            'hg19/ensembl.bed.gz',
            'hg19/ensembl.bed.gz.tbi',
            'hg19/canon_transcripts_hg19_ensembl.txt',
            'hg38/ensembl.bed.gz',
            'hg38/ensembl.bed.gz.tbi',
            'hg38/canon_transcripts_hg38_ensembl.txt',
            'canon_cancer_replacement.txt',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', script_name),
        join('scripts', 'annotate_bed.py'),
    ],
    install_requires=setup_utils.get_reqs(),
    setup_requires=['numpy'],
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)

if setup_utils.is_installing():
    print("""
--------------------------------
 {name} installation complete!
--------------------------------
Usage: {script_name} *.bam -o targqc_stats [--bed target.bed ...]'

For help in running TargQC, please see the documentation available at https://github.com/vladsaveliev/TargQC or run: targqc --help
""".format(name=name, script_name=script_name))
