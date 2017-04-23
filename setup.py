#!/usr/bin/env python
import sys
import subprocess

import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup, find_packages
from ngs_utils.setup_utils import write_version_py, run_cmdl, clean_package, find_package_files, get_reqs

name = 'TargQC'
script_name = 'targqc'
package_name = 'targqc'

if abspath(dirname(__file__)) != abspath(os.getcwd()):
    sys.stderr.write('Please, change to ' + dirname(__file__) + ' before running setup.py\n')
    sys.exit()
    

cmd = [a for a in sys.argv if not a.startswith('-')][-1]

is_installing = cmd not in ['tag', 'up', 'clean']
if is_installing:
    print('Upgrading pip and setuptools...')
    try:
        subprocess.call('pip install --upgrade pip', shell=True)
        subprocess.call('pip install --upgrade --ignore-installed setuptools', shell=True)
    except StandardError:
        sys.stderr.write('Cannot update pip and setuptools, that might cause errors '
                         'during the following intallation\n')

version = write_version_py(package_name)

if cmd == 'tag':
    run_cmdl("git tag -a %s -m \"Version %s\"" % (version, version))
    run_cmdl('git push --tags')
    sys.exit()

if cmd == 'publish':
    run_cmdl('python setup.py sdist upload')
    sys.exit()

if cmd == 'up':
    run_cmdl('git pull --recurse-submodules --rebase')
    # if first time: $ git submodule update --init --recursive
    run_cmdl('git submodule foreach "(git checkout master; git pull --rebase)"')
    sys.exit()

if cmd == 'clean':
    clean_package(package_name, '.')
    sys.exit()
    
print('Installing ' + name + ((' version ' + str(version)) if version else ''))
print('')

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
        'ngs_utils': [
        ] + find_package_files('reporting', package_name, skip_exts=['.sass', '.coffee'])
          + find_package_files('reference_data', package_name)
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
