#!/usr/bin/env python
import os
import sys
from os.path import join, isfile, abspath, dirname

from setuptools import setup, find_packages


name = 'TargQC'
script_name = 'targqc'
package_name = 'targqc'

from Utils import setup_utils
version = setup_utils.init(name, package_name, __file__)


if setup_utils.is_installing():
    bedtools_dirpath = join(setup_utils.utils_package_name, 'bedtools', 'bedtools2')
    success_compilation = setup_utils.compile_tool('BEDtools', bedtools_dirpath, [join('bin', 'bedtools')])
    if not success_compilation:
        bedtools = setup_utils.which('bedtools')
        if bedtools:
            sys.stderr.write('Compilation failed, using bedtools in $PATH: ' + bedtools + '\n')
        else:
            sys.exit(1)

    bwa_dirpath = join(package_name, 'bwa')
    success_compilation = setup_utils.compile_tool('bwa', bwa_dirpath, ['bwa'])
    if not success_compilation: sys.stderr.write('BWA has failed to compile, cannot process FastQ without BWA')


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
        'GeneAnnotation': [
            'Ensembl/biomart.tsv',
            'Ensembl/hg19/ensembl.bed.gz',
            'Ensembl/hg19/ensembl.bed.gz.tbi',
            'Ensembl/hg38/ensembl.bed.gz',
            'Ensembl/hg38/ensembl.bed.gz.tbi',
            ],
        setup_utils.utils_package_name: setup_utils.get_utils_package_files(),
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
    },
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', script_name),
        join('GeneAnnotation', 'annotate_bed.py')],
    install_requires=setup_utils.get_reqs(),
    setup_requires=[
        'numpy'
    ],
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

if sys.argv[-1] == 'install':
    print("""
--------------------------------
 {name} installation complete!
--------------------------------
Usage: {script_name} *.bam -o targqc_stats [--bed target.bed ...]'

For help in running TargQC, please see the documentation available at https://github.com/vladsaveliev/TargQC or run: targqc --help
""".format(name=name, script_name=script_name))
