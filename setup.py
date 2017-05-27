#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))

from os.path import join, isfile, abspath, dirname
from ngs_utils.setup_utils import write_version_py, find_package_files, get_reqs
from ngs_utils import setup_utils

name = 'TargQC'
script_name = 'targqc'
package_name = 'targqc'


version = setup_utils.init(package_name, package_name, __file__)


from setuptools import setup, find_packages
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
            'gender/*.bed',
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
        'bed_venn': [
            '*.html'
        ]
    },
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', script_name),
        join('scripts', 'annotate_bed.py'),
        join('scripts', 'bed_venn.py'),
    ],
    install_requires=get_reqs(),
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

print("""
--------------------------------
 {name} installation complete!
--------------------------------
Usage: {script_name} *.bam -o targqc_stats [--bed target.bed ...]'

For help in running TargQC, please see the documentation available at https://github.com/vladsaveliev/TargQC or run: targqc --help
""".format(name=name, script_name=script_name))
