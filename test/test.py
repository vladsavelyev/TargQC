import subprocess
import unittest
import os
from genericpath import getsize
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath
from datetime import datetime
from nose.plugins.attrib import attr
from collections import namedtuple

from Utils.file_utils import safe_mkdir


def call(cmdl):
    print cmdl
    subprocess.call(cmdl)

def check_call(cmdl):
    print '*' * 100
    print cmdl
    print '*' * 100
    print ''
    subprocess.check_call(cmdl)


# def run(args, fname, suf):
#     fpath = join(test_dir, fname)
#     output_fpath = join(test_dir, splitext(fname)[0] + '.anno' + (('_'+suf) if suf else '') + '.bed')
#     prev_output_fpath = None
#     if isfile(output_fpath):
#         prev_output_fpath = output_fpath + '_' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
#         os.rename(output_fpath, prev_output_fpath)
#
#     call(script + ' ' + fpath + ' ' + args + ' -o ' + output_fpath)
#
#     if prev_output_fpath:
#         os.system('diff ' + prev_output_fpath + ' ' + output_fpath)


class TestTargQC(unittest.TestCase):
    script = 'targqc'

    syn3_url = 'https://www.dropbox.com/s/lxv11oou8r316sq//syn3-chr21.tar.gz'
    bwa_url = 'https://www.dropbox.com/s/lxv11oou8r316sq/bwa.tar.gz'

    data_dir = join(dirname(__file__), 'data')

    syn3_dir = join(data_dir, 'syn3-chr21')
    bed3 = join(syn3_dir, 'NGv3.chr21.3col.bed')
    bed4 = join(syn3_dir, 'NGv3.chr21.4col.bed')
    Sample = namedtuple('Sample', 'name bam l_fastq r_fastq')
    samples = [
        Sample('syn3-normal', 'syn3-tumor.chr21.bam', 'syn3-normal_R1.fq.gz', 'syn3-normal_R2.fq.gz'),
        Sample('syn3-tumor', 'syn3-tumor.chr21.bam', 'syn3-tumor_R1.fq.gz', 'syn3-tumor_R2.fq.gz'),
    ]

    bams = [join(syn3_dir, s.bam) for s in samples]
    fastqs = [join(syn3_dir, s.l_fastq) for s in samples] + [join(syn3_dir, s.r_fastq) for s in samples]
    bwa_dir = join(data_dir, 'bwa')
    bwa_path = join(bwa_dir, 'hg19-chr21.fa')

    results_dir = safe_mkdir(join(dirname(__file__), 'results'))
    gold_standard_dir = join(join(dirname(__file__), 'gold_standard'))

    def setUp(self):
        if not exists(self.syn3_dir):
            self._download_test_data()
        if not exists(self.bwa_dir):
            self._download_bwa_data()

    def _download_test_data(self):
        os.chdir(self.data_dir)
        check_call(['wget', self.syn3_url, '-O', self.syn3_dir])
        check_call(['tar', '-xzvf', basename(self.syn3_url)])

    def _download_bwa_data(self):
        os.chdir(self.data_dir)
        check_call(['wget', self.bwa_url, '-O', self.bwa_dir])
        check_call(['tar', '-xzvf', basename(self.bwa_url)])

    def check_file(self, fpath):
        assert isfile(fpath)
        assert getsize(fpath) > 0
        if isdir(self.gold_standard_dir):
            cmp_fpath = join(self.gold_standard_dir, relpath(fpath, self.results_dir))
            if cmp_fpath:
                check_call(['diff', '-q', fpath, cmp_fpath])

    def _test(self, used_samples, output_dirname=None, bams=None, fastq=None, bed=None,
              debug=False, reannotate=False, genome=None, bwa=None, threads=None, ipython=None):
        os.chdir(self.results_dir)
        cmdl = [self.script]
        output_dir = None
        if output_dirname:
            output_dir = join(self.results_dir, output_dirname)
            cmdl.extend(['-o', output_dir])
        if bams: cmdl.extend(bams)
        if fastq: cmdl.extend(fastq)
        if bed: cmdl.extend(['--bed', bed])
        if debug: cmdl.append('--debug')
        if reannotate: cmdl.append('--reannotate')
        if genome: cmdl.extend(['-g', genome])
        if bwa: cmdl.extend(['--bwa', bwa])
        if threads: cmdl.extend(['-t', threads])
        if ipython: cmdl.extend('-s sge -q queue -r pename=smp'.split())

        check_call(cmdl)

        # Checking results
        if not output_dir:
            output_dir = join(os.getcwd(), 'targqc')
        assert isdir(output_dir)
        for fname in ['regions.tsv', 'summary.html', 'summary.tsv']:
            self.check_file(join(output_dir, fname))
            for s in used_samples:
                assert isdir(join(output_dir, s.name))
                for fname in ['regions.tsv', 'summary.html', 'summary.txt', 'summary.json']:
                    fpath = join(output_dir, s.name, fname)
                    self.check_file(fpath)
        # TODO: check line numbers and some values isntead of diff?

    def test_onesample(self):
        self._test([self.samples[0]], bams=[self.bams[0]], bed=self.bed4, genome='hg19', output_dirname='one_sample')

    def test_bed3(self):
        pass

    def test_bed4(self):
        pass

    def test_bed4_reannotate(self):
        pass

    def test_wgs(self):
        pass

    def test_ipython(self):
        pass

    def test_threads(self):
        pass

    def test_fastq(self):
        pass

    def test_nodebug(self):
        pass

    def test_hg19(self):
        pass

    def test_hg38(self):
        pass

    def test_api(self):
        pass
