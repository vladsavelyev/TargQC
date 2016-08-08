import subprocess
import unittest
import os
from genericpath import getsize
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath
from datetime import datetime
from time import sleep

from nose.plugins.attrib import attr
from collections import namedtuple
import sys

from nose import SkipTest
import nose


def info(msg=''):
    print msg
    sys.stdout.flush()

def call(cmdl):
    info(cmdl)
    subprocess.call(cmdl)

def check_call(cmdl):
    info(cmdl if isinstance(cmdl, basestring) else ' '.join(cmdl))
    subprocess.check_call(cmdl, shell=isinstance(cmdl, basestring))

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

    syn3_url = 'https://www.dropbox.com/s/lxv11oou8r316sq/TargQC_test_data/syn3-chr21.tar.gz'
    bwa_url = 'https://www.dropbox.com/s/lxv11oou8r316sq/TargQC_test_data/bwa.tar.gz'

    data_dir = join(dirname(__file__), 'data')

    syn3_dir = join(data_dir, 'syn3-chr21')
    bed3 = join(syn3_dir, 'NGv3.chr21.3col.bed')
    bed4 = join(syn3_dir, 'NGv3.chr21.4col.bed')
    Sample = namedtuple('Sample', 'name bam l_fastq r_fastq')
    samples = [
        Sample('syn3-tumor', 'syn3-tumor.bam', 'syn3-tumor_R1.fq.gz', 'syn3-tumor_R2.fq.gz'),
        Sample('syn3-normal', 'syn3-normal.bam', 'syn3-normal_R1.fq.gz', 'syn3-normal_R2.fq.gz'),
    ]

    bams = [join(syn3_dir, s.bam) for s in samples]
    fastqs = [join(syn3_dir, s.l_fastq) for s in samples] + [join(syn3_dir, s.r_fastq) for s in samples]
    bwa_dir = join(data_dir, 'bwa')
    bwa_path = join(bwa_dir, 'hg19-chr21.fa')

    results_dir = join(dirname(__file__), 'results')
    gold_standard_dir = join(join(dirname(__file__), 'gold_standard'))

    def setUp(self):
        if not exists(self.syn3_dir):
            self._download_test_data()
        if not exists(self.bwa_dir):
            self._download_bwa_data()
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _download_test_data(self):
        os.chdir(self.data_dir)
        check_call(['wget', self.syn3_url, '-O', self.syn3_dir])
        check_call(['tar', '-xzvf', basename(self.syn3_url)])

    def _download_bwa_data(self):
        os.chdir(self.data_dir)
        check_call(['wget', self.bwa_url, '-O', self.bwa_dir])
        check_call(['tar', '-xzvf', basename(self.bwa_url)])

    def check_file(self, fpath, diff_ignore_re=''):
        assert isfile(fpath)
        assert getsize(fpath) > 0
        if isdir(self.gold_standard_dir):
            cmp_fpath = join(self.gold_standard_dir, relpath(fpath, self.results_dir))
            if cmp_fpath and isfile(cmp_fpath):
                if diff_ignore_re:
                    cmdl = ['diff', '-q', '--ignore-matching-lines', diff_ignore_re, fpath, cmp_fpath]
                else:
                    cmdl = ['diff', '-q', fpath, cmp_fpath]
                check_call(cmdl)

    def _test(self, output_dirname=None, used_samples=samples, bams=None, fastq=None, bed=None,
              debug=False, reuse=False, reannotate=False, genome='hg19-chr21', bwa=None,
              threads=None, ipython=None):
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
        if reuse: cmdl.append('--reuse')
        if reannotate: cmdl.append('--reannotate')
        if genome: cmdl.extend(['-g', genome])
        if bwa: cmdl.extend(['--bwa', bwa])
        if threads: cmdl.extend(['-t', threads])
        if ipython: cmdl.extend('-s sge -q queue -r pename=smp'.split())

        info('-' * 100)
        check_call(cmdl)
        info('-' * 100)
        info('')

        self._check_results(output_dir, used_samples)

    def _check_results(self, output_dir, used_samples):
        if not output_dir:
            output_dir = join(os.getcwd(), 'targqc')
        assert isdir(output_dir)
        self.check_file(join(output_dir, 'regions.tsv'))
        self.check_file(join(output_dir, 'summary.tsv'))
        self.check_file(join(output_dir, 'summary.html'), diff_ignore_re='report_date')
        for s in used_samples:
            s_dir = join(output_dir, s.name)
            assert isdir(s_dir)
            self.check_file(join(s_dir, 'regions.tsv'))
            self.check_file(join(s_dir, 'summary.txt'))
            self.check_file(join(s_dir, 'summary.html'), diff_ignore_re='report_date')
            self.check_file(join(s_dir, 'summary.json'), diff_ignore_re='work_dir')
        # TODO: check line numbers and some values isntead of diff?

    def test_onesample(self):
        self._test('one_sample', [self.samples[0]], bams=[self.bams[0]], bed=self.bed4)

    def test_bed3(self):
        self._test('bed3', bams=self.bams, bed=self.bed3)

    def test_bed4(self):
        self._test('bed4', bams=self.bams, bed=self.bed4)

    def test_bed4_reannotate(self):
        self._test('bed4_reannotate', bams=self.bams, bed=self.bed4, reannotate=True)

    def test_wgs(self):
        self._test('wgs', bams=self.bams)

    def test_threads(self):
        self._test('threads', bams=self.bams, bed=self.bed4, threads=2)

    def test_ipython(self):
        self._test('ipython', bams=self.bams, bed=self.bed4, ipython=True, threads=2)

    def test_fastq(self):
        self._test('fastq', fastq=self.fastqs, bwa=self.bwa_path, bed=self.bed4)

    def test_debug_and_reuse(self):
        self._test('debug_and_reuse', bams=self.bams, bed=self.bed4, debug=True, reuse=True)

    def test_full_hg19(self):
        raise SkipTest

    def test_full_hg38(self):
        raise SkipTest

    # def test_api(self):
    #     import targqc
    #     import Utils.reference_data as ref
    #     from Utils.file_utils import safe_mkdir
    #     from Utils.parallel import ParallelCfg
    #
    #     genome = 'hg19-chr21'
    #     fai_fpath = ref.get_fai(genome)
    #     output_dirname = 'api'
    #     output_dir = join(self.results_dir, output_dirname)
    #     work_dir = join(output_dir, 'work')
    #     samples = [targqc.Sample(s.name,
    #          dirpath=safe_mkdir(join(output_dir, s.name)),
    #          work_dir=safe_mkdir(join(work_dir, s.name)),
    #          bam=s.bam) for s in self.samples]
    #     parallel_cfg = ParallelCfg(None, None, None, 1, None)
    #     info('-' * 100)
    #     targqc.start_targqc(work_dir, output_dir, samples, self.bed4,
    #                         parallel_cfg, self.bwa_path,
    #                         fai_fpath=fai_fpath,
    #                         genome=genome)
    #     info('-' * 100)
    #     info('')
    #
    #     info()
    #     self._check_results(output_dir, self.samples)


# nose.main()