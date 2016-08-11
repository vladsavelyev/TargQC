#!/usr/bin/env python
import os
import subprocess
import sys
from os.path import join, isfile, abspath, dirname

from setuptools import setup, find_packages
from pip.req import parse_requirements
from sys import platform as sys_platform
import platform

name = 'targqc'


try:
    from pip.download import PipSession
except ImportError:  # newer setuptools
    install_reqs = parse_requirements('requirements.txt')
else:
    install_reqs = parse_requirements('requirements.txt', session=PipSession())

reqs = [str(ir.req) for ir in install_reqs]


version = open(join(dirname(abspath(__file__)), 'VERSION.txt')).read().strip()
'''
Versioning:
1. Write each version to VERSION.txt
2. If the changes are significant, tag the release and push the new tag:
   $ python setup.py tag
'''

def write_version_py():
    version_py = os.path.join(os.path.dirname(__file__), 'targqc', 'version.py')
    try:
        import subprocess
        git_revision = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).rstrip()
    except:
        git_revision = ''
    with open(version_py, 'w') as f:
        f.write((
            '# Do not edit this file, pipeline versioning is governed by git tags\n' +
            '__version__ = \'' + version + '\'\n' +
            '__git_revision__ = \'' + git_revision + '\''))

write_version_py()

def all_required_binaries_exist(aligner_dirpath, required_binaries):
    for required_binary in required_binaries:
        if not isfile(join(aligner_dirpath, required_binary)):
            return False
    return True

def which(program):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(_fpath):
        return os.path.isfile(_fpath) and os.access(_fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def compile_tool(tool_name, dirpath, requirements):
    make_logs_basepath = join(dirpath, 'make')
    failed_compilation_flag = make_logs_basepath + '.failed'

    if not all_required_binaries_exist(dirpath, requirements):
        # making
        print 'Compiling ' + tool_name
        return_code = subprocess.call(['make', '-C', dirpath],
              stdout=open(make_logs_basepath + '.log', 'w'),
              stderr=open(make_logs_basepath + '.err', 'w'),)

        if return_code != 0 or not all_required_binaries_exist(dirpath, requirements):
            sys.stderr.write('Failed to compile ' + tool_name + ' (' + dirpath + '), '
                 'details are in ' + make_logs_basepath + '.log and make.err\n')
            open(failed_compilation_flag, 'w').close()
            return False
        try:
            os.remove(make_logs_basepath + '.log')
            os.remove(make_logs_basepath + '.err')
        except OSError as e:
            sys.stderr.write('Warning: cannot remove make_logs_basepath.* : ' + str(e) + '\n')
    return True

bedtools_dirpath = join(dirname(abspath(__file__)), 'Utils', 'bedtools', 'bedtools2')
success_compilation = compile_tool('BEDtools', bedtools_dirpath, [join('bin', 'bedtools')])
if not success_compilation:
    bedtools = which('bedtools')
    if bedtools:
        sys.stderr.write('Compilation failed, using bedtools in $PATH: ' + bedtools + '\n')
    else:
        sys.exit(1)

bwa_dirpath = join(dirname(abspath(__file__)), 'bwa')
success_compilation = compile_tool('bwa', bwa_dirpath, ['bwa'])
if not success_compilation: sys.stderr.write('BWA has failed to compile, cannot process FastQ without BWA')


def _run(_cmd):
    print('$ ' + _cmd)
    os.system(_cmd)

if sys.argv[-1] == 'publish':
    _run('python setup.py sdist upload')
    _run('python setup.py bdist_wheel upload')
    sys.exit()

if sys.argv[-1] == 'tag':
    _run('git tag -a %s -m "Version %s"' % (version, version))
    _run('git push --tags')
    sys.exit()

if sys.argv[-1] == 'up':
    _run('cd Utils && git up && cd ..')
    _run('cd GeneAnnotation && git up && cd ..')
    _run('cd MultiQC && git up && cd ..')
    _run('git up')
    sys.exit()

if sys.argv[-1] == 'install':
    print('''-----------------------------------
 Installing TargQC version {}
-----------------------------------
'''.format(version))

# sambamba = join('Utils', 'sambamba', 'build', 'sambamba')
# if not isfile(sambamba):
#     _run('cd Utils/sambamba; make sambamba-ldmd2-64; cd ../..')
#     if not isfile(sambamba):
#         sys.stderr.write('Error: could not compile sambamba, exiting.\n')
#         sys.exit(1)

def sambamba_executable():
    sambamba_dirpath = join('Utils', 'sambamba_binaries')
    if 'darwin' in sys_platform:
        path = join(sambamba_dirpath, 'sambamba_osx')
    elif 'redhat' in platform.dist():
        path = join(sambamba_dirpath, 'sambamba_centos')
    else:
        path = join(sambamba_dirpath, 'sambamba_lnx')
    if isfile(path):
        return path
    elif isfile(path + '.gz'):
        print('gunzipping sambamba ' + path + '.gz')
        os.system('gunzip ' + path + '.gz')
        return path
    else:
        sys.stderr.write('Error: could not find sambamba ' + path + '(.gz)')

_run('cd MultiQC && python setup.py install && cd ..')

setup(
    name=name,
    version=version,
    author='Vlad Saveliev and Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='Target coverage evaluation tool',
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
        'Utils': [
            'reference_data/fai/*.fai',
            'reporting/static/*.js',
            'reporting/static/*/*.js',
            'reporting/static/*.css',
            'reporting/static/*/*.css',
            'reporting/static/*.json',
            'reporting/static/*/*.json',
            'reporting/static/*.png',
            'reporting/static/*/*.png',
            'reporting/static/*.pxm',
            'reporting/static/*/*.pxm',
            'reporting/*.html',
            'reporting/*.json',
            os.path.relpath(sambamba_executable(), 'Utils'),
            'bedtools/bedtools2/bin/*',
            'tools/*.sh',
        ],
        'targqc': [
            'bedops/bedops_*',
            'qualimap/*/qualimap',
            'qualimap/*/qualimap.jar',
            'qualimap/*/lib/*.jar',
            'qualimap/*/scripts/*.jar',
            'qualimap/*/species/*.jar',
            'picard/picard/*.jar',
            'gender/*.bed',
        ],
        'bwa': [
            'bwa',
        ]
    },
    include_package_data=True,
    zip_safe=False,
    scripts=['scripts/' + name, 'GeneAnnotation/annotate_bed.py'],
    install_requires=reqs,
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

if sys.argv[-1] == 'install':
    print("""
--------------------------------
 TargQC installation complete!
--------------------------------
Usage: {name} *.bam -o targqc_stats [--bed target.bed ...]'

For help in running TargQC, please see the documentation available at https://github.com/vladsaveliev/TargQC or run: targqc --help
""".format(name=name))
