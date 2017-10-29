#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))

from os.path import join, isfile, abspath, dirname

name = 'TargQC'
script_name = 'targqc'
package_name = 'targqc'


from setuptools import setup, find_packages

try:
    from utilz import setup_utils
except:
    version = open('VERSION.txt').read().strip().split('\n')[0]
    setup(version=version)  # For conda-build jinja context to read the version
else:
    version = setup_utils.init(package_name, package_name, __file__)
    from utilz.setup_utils import write_version_py, find_package_files, get_reqs
    from utilz import setup_utils
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
            ]
        },
        include_package_data=True,
        zip_safe=False,
        scripts=[
            join('scripts', script_name),
            join('scripts', 'annotate_bed.py'),
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