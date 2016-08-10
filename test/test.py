import subprocess
import unittest
import os
from genericpath import getsize, getmtime
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath
from datetime import datetime
from time import sleep

import shutil
from nose.plugins.attrib import attr
from collections import namedtuple
import sys

from nose import SkipTest
import nose


def info(msg=''):
    sys.stdout.write(msg + '\n')

def call(cmdl):
    info(cmdl)
    subprocess.call(cmdl)

def check_call(cmdl):
    info(cmdl if isinstance(cmdl, basestring) else ' '.join(cmdl))
    subprocess.check_call(cmdl, shell=isinstance(cmdl, basestring))


class TestTargQC(unittest.TestCase):
    script = 'targqc'

    syn3_url = 'http://quast.bioinf.spbau.ru/static/chr21.tar.gz'
    bwa_url = 'http://quast.bioinf.spbau.ru/static/bwa.tar.gz'

    data_dir = join(dirname(__file__), 'data')

    syn3_dir = join(data_dir, 'chr21')
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

    remove_work_dir_on_success = False

    def setUp(self):
        if not isdir(self.data_dir):
            os.mkdir(self.data_dir)
        if not isdir(self.syn3_dir):
            info(self.syn3_dir + ' does not exist, downloading test data')
            cur_dir = os.getcwd()
            os.chdir(self.data_dir)
            check_call(['wget', self.syn3_url])
            check_call(['tar', '-xzvf', basename(self.syn3_url)])
            os.chdir(cur_dir)
        if not isdir(self.bwa_dir):
            info(self.bwa_dir + ' does not exist, downloading bwa reference data')
            cur_dir = os.getcwd()
            os.chdir(self.data_dir)
            check_call(['wget', self.bwa_url])
            check_call(['tar', '-xzvf', basename(self.bwa_url)])
            os.chdir(cur_dir)
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _check_file(self, fpath, diff_ignore_re=''):
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
              debug=False, reuse_intermediate=False, reuse_output_dir=False, reannotate=False,
              genome='hg19-chr21', bwa=None, threads=None, ipython=None, keep_work_dir=True):
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
        if reuse_intermediate:
            reuse_output_dir = True
            cmdl.append('--reuse')
        if reannotate: cmdl.append('--reannotate')
        if genome: cmdl.extend(['-g', genome])
        if bwa: cmdl.extend(['--bwa', bwa])
        if threads: cmdl.extend(['-t', str(threads)])
        if ipython: cmdl.extend('-s sge -q queue -r pename=smp -r --local'.split())
        if keep_work_dir: cmdl.append('--keep-work-dir')

        output_dir = output_dir or self._default_output_dir()
        if exists(output_dir) and reuse_output_dir is False:
            last_changed = datetime.fromtimestamp(getmtime(output_dir))
            prev_output_dir = output_dir + '_' + last_changed.strftime("%Y_%m_%d_%H_%M_%S")
            os.rename(output_dir, prev_output_dir)

        info('-' * 100)
        check_call(cmdl)
        info('-' * 100)
        info('')

        self._check_results(output_dir, used_samples)

        if self.remove_work_dir_on_success and not reuse_intermediate and not reuse_output_dir:
            work_dir = join(output_dir, 'work')
            if not isdir(work_dir):
                info('Work dir for run ' + output_dirname + ' does not exist under ' + work_dir)
            else:
                shutil.rmtree(work_dir)

    def _default_output_dir(self):
        return join(os.getcwd(), 'targqc')

    def _check_results(self, output_dir, used_samples):
        assert isdir(output_dir)
        self._check_file(join(output_dir, 'regions.tsv'))
        self._check_file(join(output_dir, 'summary.tsv'))
        self._check_file(join(output_dir, 'summary.html'), diff_ignore_re='report_date')
        for s in used_samples:
            s_dir = join(output_dir, s.name)
            assert isdir(s_dir)
            self._check_file(join(s_dir, 'regions.tsv'))
            self._check_file(join(s_dir, 'summary.txt'))
            self._check_file(join(s_dir, 'summary.html'), diff_ignore_re='report_date')
            self._check_file(join(s_dir, 'summary.json'), diff_ignore_re='work_dir')
        # TODO: check line numbers and some values isntead of diff?

    def test_01_simple(self):
        self._test('simple', [self.samples[0]], bams=[self.bams[0]], bed=self.bed4)

    def test_02_bed3(self):
        self._test('bed3', bams=self.bams, bed=self.bed3)

    def test_03_bed4(self):
        self._test('bed4', bams=self.bams, bed=self.bed4)

    def test_04_bed4_reannotate(self):
        self._test('bed4_reannotate', bams=self.bams, bed=self.bed4, reannotate=True)

    def test_05_wgs(self):
        self._test('wgs', bams=self.bams)

    def test_06_threads(self):
        self._test('threads', bams=self.bams, bed=self.bed4, threads='2')

    # def test_07_ipython(self):
    #     self._test('ipython', bams=self.bams, bed=self.bed4, ipython=True, threads='2')

    def test_08_fastq(self):
        self._test('fastq', fastq=self.fastqs, bwa=self.bwa_path, bed=self.bed4)

    def test_09_debug_and_reuse(self):
        self._test('debug_and_reuse', bams=self.bams, bed=self.bed4, debug=True, reuse_intermediate=True)

    def test_09_reuse_output(self):
        self._test('reuse_output_dir', bams=self.bams, bed=self.bed4, reuse_output_dir=True)

    def test_10_full_hg19(self):
        raise SkipTest

    def test_11_full_hg38(self):
        raise SkipTest

    def test_12_api(self):
        import targqc
        import Utils.reference_data as ref
        from Utils.file_utils import safe_mkdir
        from Utils.parallel import ParallelCfg

        genome = 'hg19-chr21'
        fai_fpath = ref.get_fai(genome)
        output_dirname = 'api'
        output_dir = join(self.results_dir, output_dirname)
        work_dir = join(output_dir, 'work')
        samples = sorted([targqc.Sample(s.name,
             dirpath=safe_mkdir(join(output_dir, s.name)),
             work_dir=safe_mkdir(join(work_dir, s.name)),
             bam=join(self.syn3_dir, s.bam))
                   for s in self.samples], key=lambda _s: _s.key_to_sort())
        parallel_cfg = ParallelCfg(None, None, None, 1, None)
        info('-' * 100)
        targqc.start_targqc(work_dir, output_dir, samples, self.bed4,
                            parallel_cfg, self.bwa_path,
                            fai_fpath=fai_fpath,
                            genome=genome)
        info('-' * 100)
        info('')

        info()
        self._check_results(output_dir, self.samples)


# nose.main()