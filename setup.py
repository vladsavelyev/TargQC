#!/usr/bin/env python

from setuptools import setup, find_packages

version = '0.1'

print("""-----------------------------------
 Installing TargQC version {}
-----------------------------------
""".format(version))

setup(
    name='targqc',
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Target coverage evaluation tool',
    long_description=__doc__,
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/TargQC',
    download_url='https://github.com/vladsaveliev/TargQC/releases',
    license='GPLv3',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    scripts=['scripts/targqc', 'scripts/targqc_multi'],
    install_requires=[
        'pybedtools',
        'pysam',
        'ipython-cluster-helper',
        # 'pyyaml',
        # 'click',
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
)

print("""
--------------------------------
 TargQC installation complete!
--------------------------------
For help in running TargQC, please see the documentation available
at https://github.com/vladsaveliev/TargQC or run: targqc --help
""")
